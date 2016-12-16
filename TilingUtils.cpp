#include "TilingUtils.h"

#include <iostream>
#include <set>
#include <numeric> // iota

#include <igl/components.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>
#include "offsetSurface.h"
#include "glob_defs.h"


extern "C" 
{
#include <tourtre.h>
}

using namespace std;

namespace {

int nComponents(const Eigen::VectorXi &comps) {
  std::set<int> u;
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
  void getNeighbors(size_t idx, std::vector<size_t> &n) const {
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
  std::vector<size_t> nbrs_buf;

  tdat->getNeighbors(idx, nbrs_buf);
  for (int i = 0; i < nbrs_buf.size(); ++i) {
    nbrs[i] = nbrs_buf[i];
  }
  return nbrs_buf.size();
}

void extractTreeVerts(ctBranch *b, std::vector<TilingUtils::ConnectedComponent> &ccs, tourtre_data *tdat) {
  // At each node in the tree, generate the offset surface.
  int off_vert = b->extremum;
  double off_value = tdat->_H[off_vert];

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi O;

  /*
  offsetSurface::generateOffsetSurface_naive(tdat->_V, tdat->_T, tdat->_F, tdat->_O,
                                             off_value, V, F, O);
                                             */
}

void outputTree(ctBranch* b, tourtre_data* tdat, std::vector<int> &nodes,
                std::set<double> &unique_offsets) {
  if ( (tdat->_H(b->saddle) > 0 && tdat->_H(b->saddle) < 1)) {
    unique_offsets.insert(tdat->_H(b->extremum));
    unique_offsets.insert(tdat->_H(b->saddle));
    printf("(%d:%f %d:%f", b->extremum, tdat->_H(b->extremum),
           b->saddle, tdat->_H(b->saddle));
  }
  nodes.push_back(b->extremum);
  nodes.push_back(b->saddle);

	for ( ctBranch * c = b->children.head; c != NULL; c = c->nextChild ){
    printf(" ");
		outputTree( c, tdat, nodes, unique_offsets);
	}
	
	printf(")");
}

} // namespace

namespace TilingUtils {

ConnectedComponent::ConnectedComponent(double offset,
                                       const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
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
  this->O.resize(vertices_used.size(), 3);
  this->F.resize(faces_used.size(), 3);

  // Create mapping and adjust vertices.
  map<int, int> original_to_new_vertex;
  int new_vertex = 0;
  for (int original_vertex : vertices_used) {
    this->V.row(new_vertex) = V.row(original_vertex);
    this->O.row(new_vertex) = O.row(original_vertex);
    original_to_new_vertex[original_vertex] = new_vertex++;
  }

  // Adjust faces.
  int face_index = 0;
  for (int i : faces_used) {
    for (int j = 0; j < F.cols(); j++)
      this->F(face_index, j) = original_to_new_vertex[F(i, j)];
    face_index++;
  }
}

void ConnectedComponent::render() const {
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(this->V, this->F);
  viewer.launch();
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
    components.push_back(ConnectedComponent(offset, V, F, vertices_used));
  }

  return components;
}

void addVerticesToContour(int vertex_index,
                          set<int>& visited, set<int> &vertices_used,
                          const Eigen::MatrixXi &F, const Eigen::VectorXi &O) {
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
      else if (O(neighbor_vertex) != GLOBAL::original_marker)
        continue;

      visited.insert(neighbor_vertex);
      vertices_used.insert(neighbor_vertex);
      addVerticesToContour(neighbor_vertex, visited, vertices_used, F, O);
    }
  }
}

vector<set<int> > getContourVertices(const Eigen::MatrixXi &F,
                                     const Eigen::VectorXi &O) {
  vector<set<int> > all_contours;

  set<int> visited;
  for (int i = 0; i < O.rows(); ++i) {
    if (O(i) != GLOBAL::original_marker)
      continue;
    else if (visited.find(i) != visited.end())
      continue;

    visited.insert(i);

    all_contours.push_back(set<int> ());

    set<int> &contour_vertices = all_contours.back();
    contour_vertices.insert(i);
    addVerticesToContour(i, visited, contour_vertices, F, O);
  }

  return all_contours;
}

set<int> getContoursUsed(const ConnectedComponent &component,
                         vector<set<int> > &contours,
                         const Eigen::MatrixXd &V) {
  set<int> used;

  for (int i = 0; i < component.V.rows(); i++) {

    for (int j = 0; j < contours.size(); j++) {
      set<int> &contour = contours[j];

      for (int vertex : contour) {
        double dist = 0.0;
        for (int k = 0; k < V.cols(); k++) {
          double tmp = component.V(i, k) - V(vertex, k);
          dist += tmp * tmp;
        }
        if (dist <= 1e-3)
          used.insert(j);
      }
    }
  }

  cout << endl;
  for (int a : used)
    cout << a << " ";
  cout << endl;

  return used;
}

std::vector<ConnectedComponent> allPossibleTiles(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO, const Eigen::VectorXd &H) {
  // Use libtourtre to find all possible saddle points, which leads to all possible
  // connected components.
  // Create data struct
  tourtre_data tdat(TV, TF, TT, H, TO);
  // Create sorted indices.
  std::vector<size_t> orders;
  orders.resize(TV.rows());
  std::iota(orders.begin(), orders.end(), 0);
  std::sort(orders.begin(), orders.end(), 
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

  std::vector<int> nodes;
  std::set<double> unique_offsets;
  outputTree(root, &tdat, nodes, unique_offsets);
  Eigen::MatrixXd pts, pts_cols;
  Eigen::VectorXd pts_cols_i;
  pts.resize(nodes.size(), 3);
  pts_cols_i.resize(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    pts.row(i) = TV.row(nodes[i]);
    pts_cols_i(i) = H(nodes[i]);
  }
  igl::jet(pts_cols_i, true, pts_cols);

  Eigen::MatrixXd cols;
  igl::jet(TO, true, cols);

  printf("TV: %lu TF:%lu pts:%lu\n", TV.rows(), TF.rows(), pts.rows());


  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(TV, TF);
  viewer.data.set_colors(cols);
  viewer.data.set_points(pts, pts_cols);
  bool display_mesh = true;
  viewer.callback_key_down = [&](igl::viewer::Viewer& v, unsigned char key, int modifier) {
      if (key == ' ') {
        display_mesh = !display_mesh;
        if (display_mesh) {
          v.data.set_mesh(TV, TF);
          v.data.set_colors(cols);
          v.data.add_points(pts, pts_cols);
        } else {
          v.data.clear();
          v.data.add_points(pts, pts_cols);
        }
        return true;
      }
      return false;
    };
  viewer.launch();

  vector<ConnectedComponent> all_ccs;
  vector<set<int> > cc_contours;

  vector<set<int> > contours = getContourVertices(TF, TO);

  for (double offset: unique_offsets) {
    if (offset < 0.1)
      continue;
    cout << offset << endl;

    Eigen::MatrixXd offsetV;
    Eigen::MatrixXi offsetF;
    Eigen::VectorXi offsetO;
    OffsetSurface::generateOffsetSurface_naive(TV, TT, TO, H,
                                               offset,
                                               offsetV, offsetF, offsetO);
    viewer.data.clear();
    viewer.data.set_mesh(offsetV, offsetF);
    viewer.launch();

    vector<ConnectedComponent> components = getConnectedComponents(offsetV,
                                                                   offsetF,
                                                                   offsetO,
                                                                   offset);
    for (const ConnectedComponent &component : components) {
      bool is_unique = true;

      set<int> tmp = getContoursUsed(component, contours, TV);
      for (set<int> &cc_contour : cc_contours)
        is_unique &= (tmp != cc_contour);

      if (is_unique) {
        cc_contours.push_back(tmp);
        all_ccs.push_back(component);
        component.render();
      }
    }
  }

  /*
  int num_cc = -1;
  igl::components(offsetF, comp);
  while(num_cc != 1) {
  }
  */
  return all_ccs;
}

} // namespace TilingUtils
