#include "TilingUtils.h"

#include <iostream>
#include <set>
#include <numeric> // iota

#include <igl/components.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>
#include "offsetSurface.h"
#include "glob_defs.h"
#include "Helpers.h"
#include "viewTetMesh.h"

extern "C"
{
#include <tourtre.h>
}

using namespace std;

namespace {

int nComponents(const Eigen::VectorXi &comps) {
  set<int> u;
  for (int i = 0; i < comps.rows(); ++i) {
    u.insert(comps[i]);
  }
  return u.size();
}

struct tourtre_data {
  const Eigen::MatrixXd &_V;
  const Eigen::MatrixXi &_F;
  const Eigen::MatrixXi &_T;
  const Eigen::VectorXd &_H;
  const Eigen::VectorXi &_O;

  tourtre_data(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
              const Eigen::MatrixXi &T, const Eigen::VectorXd &H,
              const Eigen::VectorXi &O)
      : _V(V), _F(F), _T(T), _H(H), _O(O) {}

  // Needs function to get neighbors.
  void getNeighbors(size_t idx, vector<size_t> &n) const {
    for (int i = 0; i < _T.rows(); ++i) {
      bool found = false;
      for (int j = 0; j < _T.cols(); ++j) {
        if (_T(i, j) == idx) {
          found = true;
          break;
        }
      }

      // This is a tet that includes idx
      if (found) {
        // Add all the vertices in this tet that are not itself
        for (int j = 0; j < _T.cols(); ++j) {
          if (_T(i, j) != idx) {
            // SLOW (but okay): see if we've already added it (no dups)
            bool already_included = false;
            for (int k = 0; k < n.size(); ++k) {
              if (n[i] == _T(i, j)) {
                already_included = true;
                break;
              }
            }
            // Finally, add it.
            if (!already_included) {
              n.push_back(_T(i, j));
            }
          }
        }
      }
    }
  }
}; // struct tourtre_data

// Helper functions for tourtre_data
//
// Gets the value of given idx
double value (size_t idx, void * d) {
  tourtre_data* tdat = static_cast<tourtre_data*>(d);
  return tdat->_H(idx);
}

// Gets the neighbors of idx. Input pointer nbrs has already been initialized
// with at least 256 elements (change with ct_maxValence). Returns the number
// of valid entries.
size_t neighbors(size_t idx, size_t *nbrs, void* d) {
  tourtre_data* tdat = static_cast<tourtre_data*>(d);
  vector<size_t> nbrs_buf;

  tdat->getNeighbors(idx, nbrs_buf);
  for (int i = 0; i < nbrs_buf.size(); ++i) {
    nbrs[i] = nbrs_buf[i];
  }
  return nbrs_buf.size();
}

void outputTree(ctBranch* b, tourtre_data* tdat, vector<int> &nodes,
                set<double> &unique_offsets) {
  if ( (tdat->_H(b->saddle) > 0 && tdat->_H(b->saddle) < 1)) {
    unique_offsets.insert(tdat->_H(b->extremum));
    unique_offsets.insert(tdat->_H(b->saddle));
    // printf("(%d:%f %d:%f)\n", b->extremum, tdat->_H(b->extremum),
    //        b->saddle, tdat->_H(b->saddle));
  }
  nodes.push_back(b->extremum);
  nodes.push_back(b->saddle);

	for ( ctBranch * c = b->children.head; c != NULL; c = c->nextChild ){
		outputTree( c, tdat, nodes, unique_offsets);
	}
}

} // namespace

namespace TilingUtils {

ConnectedComponent::ConnectedComponent(double offset,
                                       const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
                                       const Eigen::VectorXi &M,
                                       const set<int> &vertices_used) :
    offsetVal(offset) {
  // Get all the faces used by the vertices.
  set<int> faces_used;
  for (int i = 0; i < F.rows(); i++) {
    bool use = true;
    for (int j = 0; j < F.cols(); j++)
      use &= (vertices_used.find(F(i, j)) != vertices_used.end());
    if (use)
      faces_used.insert(i);
  }

  // Fix to appropriate sizes.
  this->V.resize(vertices_used.size(), 3);
  this->F.resize(faces_used.size(), 3);
  this->M.resize(vertices_used.size());

  // Create mapping and adjust vertices.
  map<int, int> original_to_new_vertex;
  int new_vertex = 0;
  for (int original_vertex : vertices_used) {
    this->V.row(new_vertex) = V.row(original_vertex);
    this->M(new_vertex) = M(original_vertex);
    original_to_new_vertex[original_vertex] = new_vertex++;
  }

  // Adjust faces.
  int face_index = 0;
  for (int i : faces_used) {
    for (int j = 0; j < F.cols(); j++)
      this->F(face_index, j) = original_to_new_vertex[F(i, j)];
    face_index++;
  }

  // This was generated by the following contours.
  for (int i = 0; i < this->M.rows(); i++) {
    if (this->M(i) > GLOBAL::original_marker)
      this->used.insert(this->M(i));
  }
}

void dfs(int vertex_index, set<int> &visited, set<int> &vertices_used,
         const Eigen::MatrixXi &F) {
  for (int i = 0; i < F.rows(); i++) {
    bool explore = false;
    for (int j = 0; j < F.cols(); j++)
      explore |= (F(i, j) == vertex_index);
    if (!explore)
      continue;

    for (int j = 0; j < F.cols(); j++) {
      int neighbor_vertex = F(i, j);
      if (visited.find(neighbor_vertex) != visited.end())
        continue;

      visited.insert(neighbor_vertex);
      vertices_used.insert(neighbor_vertex);
      dfs(neighbor_vertex, visited, vertices_used, F);
    }
  }
}

vector<ConnectedComponent> getConnectedComponents(const Eigen::MatrixXd &V,
                                                  const Eigen::MatrixXi &F,
                                                  const Eigen::VectorXi &O,
                                                  double offset) {
  vector<ConnectedComponent> components;

  set<int> visited;
  for (int i = 0; i < V.rows(); i++) {
    if (visited.find(i) != visited.end())
      continue;
    visited.insert(i);

    set<int> vertices_used;
    dfs(i, visited, vertices_used, F);
    components.push_back(ConnectedComponent(offset, V, F, O, vertices_used));
  }

  return components;
}

vector<ConnectedComponent> allPossibleTiles(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO, const Eigen::VectorXd &H) {
  // Use libtourtre to find all possible saddle points, which leads to all possible
  // connected components.
  // Create data struct
  tourtre_data tdat(TV, TF, TT, H, TO);
  // Create sorted indices.
  vector<size_t> orders;
  orders.resize(TV.rows());
  iota(orders.begin(), orders.end(), 0);
  sort(orders.begin(), orders.end(), 
            [H](size_t i1, size_t i2) { return H(i1) < H(i2); });
  // Init libtourtre.
  ctContext *ctx = ct_init(
      TV.rows(), // # vertices
      &(orders.front()), // address of the front of the vector gives same as a C array
      &value, // function that returns the value of a given index
      &neighbors, // function that returns the number of & neighbors of a given index
      &tdat // data for callbacks
      );

  // Extract information.
  ct_sweepAndMerge(ctx);
  ctBranch* root = ct_decompose(ctx);
  ct_cleanup(ctx);

  vector<int> nodes;
  set<double> unique_offsets;
  outputTree(root, &tdat, nodes, unique_offsets);
  Eigen::MatrixXd pts, pts_cols;
  Eigen::VectorXd pts_cols_i;
  pts.resize(nodes.size(), 3);
  pts_cols_i.resize(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    pts.row(i) = TV.row(nodes[i]);
    pts_cols_i(i) = H(nodes[i]);
  }

  vector<ConnectedComponent> unique_components;

  // Go through all offsets and generate surfaces.
  for (double offset: unique_offsets) {
    // TODO: fix this.
    // There is a small offset that creates a crappy surface.
    if (offset < 0.1)
      continue;

    Eigen::MatrixXd offsetV;
    Eigen::MatrixXi offsetF;
    Eigen::VectorXi offsetO;

    OffsetSurface::generateOffsetSurface_naive(TV, TT, TO, H,
                                               offset,
                                               offsetV, offsetF, offsetO);

    vector<ConnectedComponent> components = getConnectedComponents(offsetV,
                                                                   offsetF,
                                                                   offsetO,
                                                                   offset);

    for (ConnectedComponent &component : components) {
      bool is_unique = true;
      for (const ConnectedComponent& unique_component : unique_components)
        is_unique &= (component.used != unique_component.used);

      if (is_unique)
        unique_components.push_back(component);
    }
  }

  return unique_components;
}

} // namespace TilingUtils