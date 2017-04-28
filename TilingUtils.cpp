#include "TilingUtils.h"

#include <iostream>
#include <set>
#include <numeric> // iota

#include <igl/components.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>

#include "Helpers.h"
#include "curvatureFlow.h"
#include "glob_defs.h"
#include "offsetSurface.h"
#include "viewTetMesh.h"

using namespace std;

namespace TilingUtils {

ConnectedComponent::ConnectedComponent(double offset,
                                       const Eigen::MatrixXd &V,
                                       const Eigen::MatrixXi &F,
                                       const Eigen::VectorXi &M,
                                       const set<int> &vertices_used) :
    offsetVal(offset) {
  int numDropped = 0;
  // Get all the faces used by the vertices.
  set<int> faces_used;
  for (int i = 0; i < F.rows(); i++) {
    bool use = true;
    for (int j = 0; j < F.cols(); j++)
      use &= (vertices_used.find(F(i, j)) != vertices_used.end());
    if (use)
      faces_used.insert(i);
    else
      numDropped++;
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
      this->contours_used.insert(this->M(i));
  }
}

void dfs(int vertex_index, set<int> &visited, set<int> &vertices_used,
         const map<int, vector<int> > &graph) {
  vertices_used.insert(vertex_index);
  for (int neighbor_vertex : graph.at(vertex_index)) {
    if (visited.find(neighbor_vertex) != visited.end())
      continue;

    visited.insert(neighbor_vertex);
    dfs(neighbor_vertex, visited, vertices_used, graph);
  }
}

map<int, vector<int> > getGraphFromFaces(const Eigen::MatrixXi &F) {
  map<int, vector<int> > graph;

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = j+1; k < 3; k++) {
        int u = F(i, j);
        int v = F(i, k);
        graph[u].push_back(v);
        graph[v].push_back(u);
      }
    }
  }

  return graph;
}

vector<ConnectedComponent> getConnectedComponents(const Eigen::MatrixXd &V,
                                                  const Eigen::MatrixXi &F,
                                                  const Eigen::VectorXi &O,
                                                  double offset) {
  vector<ConnectedComponent> components;
  map<int, vector<int> > graph = getGraphFromFaces(F);

  set<int> visited;
  for (int i = 0; i < V.rows(); i++) {
    if (visited.find(i) != visited.end())
      continue;
    visited.insert(i);

    set<int> vertices_used;
    dfs(i, visited, vertices_used, graph);
    components.push_back(ConnectedComponent(offset, V, F, O, vertices_used));
  }

  return components;
}

vector<ConnectedComponent> allPossibleTiles(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO, const Eigen::VectorXd &H) {

  vector<ConnectedComponent> unique_components;

  // Go through all offsets and generate surfaces.
  // for (double offset: unique_offsets) {
  // for (double offset = 0.9; offset >= 0.1; offset -= 0.1) {
  for (double offset = 0.1; offset <= 1.0; offset += 0.1) {
    Eigen::MatrixXd offsetV;
    Eigen::MatrixXi offsetF;
    Eigen::VectorXi offsetO;

    OffsetSurface::generateOffsetSurface_naive(TV, TT, TO, H,
                                               offset, offsetV, offsetF, offsetO);

    vector<ConnectedComponent> components = getConnectedComponents(offsetV,
                                                                   offsetF,
                                                                   offsetO,
                                                                   offset);

    // Many components have similar topologies - pick the best.
    for (ConnectedComponent &component : components) {
      int index = -1;

      for (int i = 0; i < unique_components.size(); i++) {
        if (unique_components[i].contours_used == component.contours_used)
          index = i;
      }

      // If this is a new combination, used it.
      if (index == -1)
        unique_components.push_back(component);
    }
  }

  return unique_components;
}

map<set<int>, ConnectedComponent> possibleTileMap(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF,
    const Eigen::MatrixXi &TT, const Eigen::VectorXi &TO,
    const Eigen::VectorXd &H) {
  map<set<int>, ConnectedComponent> result;

  // Go through all offsets and generate surfaces. Start with 1 and go to
  // lower so we emphasize larger connected components.
  //for (double offset = 0.1; offset <= 1.0; offset += 0.1) {
  for (double offset = 0.9; offset >= 0.1; offset -= 0.1) {
    Eigen::MatrixXd offsetV;
    Eigen::MatrixXi offsetF;
    Eigen::VectorXi offsetO;

    // OffsetSurface::generateOffsetSurface_naive(TV, TT, TO, H,
    //                                            offset, offsetV, offsetF, offsetO);

    OffsetSurface::marchingOffsetSurface(TV, TT, TO, H,
                                         offset, offsetV, offsetF, offsetO);

    vector<ConnectedComponent> components =
        getConnectedComponents(offsetV, offsetF, offsetO, offset);

    //printf("At offset %lf, here's the comp\n", offset);
    // Many components have similar topologies - pick the ones not found yet.
    for (ConnectedComponent &component : components) {
      /*
      printf("Used: ");
      for (auto s : component.contours_used) {
        printf("%d ", s);
      }
      printf("\n");
      */
      if (result.find(component.contours_used) == result.end()) {
        result[component.contours_used] = component;
        //Helpers::viewTriMesh(offsetV, offsetF, offsetO);
      }
    }
  }

  // Maybe haven't created a single component. If so, get larger until we do.
  double offset = 0.9;
  while (result.size() < 1) {
    // Get closer to 1 without ever getting there.
    offset = (1 + offset) / 2;

    Eigen::MatrixXd offsetV;
    Eigen::MatrixXi offsetF;
    Eigen::VectorXi offsetO;

    // OffsetSurface::generateOffsetSurface_naive(TV, TT, TO, H,
    //                                            offset, offsetV, offsetF, offsetO);

    OffsetSurface::marchingOffsetSurface(TV, TT, TO, H,
                                         offset, offsetV, offsetF, offsetO);

    vector<ConnectedComponent> components = 
        getConnectedComponents(offsetV, offsetF, offsetO, offset);
    for (ConnectedComponent &component : components) {
      if (result.find(component.contours_used) == result.end()) {
        result[component.contours_used] = component;
      }
    }
  }

  return result;
}

} // namespace TilingUtils
