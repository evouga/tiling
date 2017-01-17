#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <Eigen/Core>
#include <igl/writeOFF.h>

#include "Tiler.h"
#include "SliceStack.h"
#include "glob_defs.h"
#include "viewTetMesh.h"
#include "TilingUtils.h"
#include "Helpers.h"
#include "curvatureFlow.h"

using namespace std;
using namespace Tiler;
using namespace TilingUtils;

// Holds possible configurations from index 1 and on.
vector<Tile*> generated[1005];

// Generates contour ids.
int total = 0;

// Bot.
Eigen::MatrixXd botV;
Eigen::MatrixXi botF;
Eigen::VectorXi botM;

// Top.
Eigen::MatrixXd topV;
Eigen::MatrixXi topF;
Eigen::VectorXi topM;

// Tets.
Eigen::MatrixXd tetV;
Eigen::MatrixXi tetT;
Eigen::MatrixXi tetF;
Eigen::VectorXi tetM;

// Heat.
Eigen::VectorXd H;

// The tile falls into the same "bin" if they have the same parent connectivity.
string terribleHashFunction (const Tile *tile) {
  vector<set<int> > connected = tile->getUpperConnected();
  sort(connected.begin(), connected.end());

  stringstream ss;
  for (const set<int> &component : connected) {
    for (int node : component)
      ss << node << " ";
    ss << endl;
  }
  return ss.str();
}

double energy (Tile *tile) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi O;
  getTileMesh(tile, V, F, O, false);

  // Genus computation.
  if ((V.rows() - 0.5 * F.rows() - 2) / 2 != 0)
    return INFINITY;

  Eigen::MatrixXd unused;
  double score = biharmonic(V, F, O, unused);

  if (isnan(score))
    return INFINITY;

  return score;
}

void viewTile(Tile *tile) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi M;
  getTileMesh(tile, V, F, M, true);
  Helpers::viewTriMesh(V, F, M);
}

// NOTE: the tet-mesh should be created before calling this function.
ConnectedComponent componentFromContours(SliceStack &ss,
                                         const set<int> &contours_allowed) {
  ss.computeLaplace(tetV, tetT, tetF, tetM, H, contours_allowed);

  vector<ConnectedComponent> ccs = allPossibleTiles(tetV, tetF, tetT, tetM, H);

  if (ccs.size() == 0)
    cout << "not found" << endl;

  ConnectedComponent result = ccs[0];
  double result_score = -1.0;

  // TODO(bradyz): find a scoring mechanism.
  for (const ConnectedComponent &cc : ccs) {
    if (cc.contours_used == contours_allowed) {
      result = cc;
      result_score = 1.0;
    }
  }

  // TODO: throw a more detailed exception.
  if (result_score < 0.0) {
    cout << "not found" << endl;
    return ccs[1000];
  }

  return result;
}

// The map should contain a mapping from contours used to the components
// that have this topology.
void populateComponentMeshes(const map<set<int>, vector<Component*> > &mapping,
                             SliceStack &ss) {
  cout << "Generating meshes for " << mapping.size() << " components." << endl;

  for (auto it : mapping) {
    ConnectedComponent tmp = componentFromContours(ss, it.first);

    for (Component *component : it.second) {
      component->V = tmp.V;
      component->F = tmp.F;
      component->M = tmp.M;
    }
  }
}

vector<set<int> > getHeatFlowValidComponents(SliceStack &ss, int level) {
  ss.triangulateSlice(level, GLOBAL::TRI_AREA,
                      botV, botF, topV, topF, botM, topM);
  ss.tetrahedralizeSlice(botV, botF, topV, topF,
                         botM, topM, tetV, tetT, tetF, tetM);
  ss.computeLaplace(tetV, tetT, tetF, tetM, H);

  vector<set<int> > result;

  for (const ConnectedComponent &cc : allPossibleTiles(tetV, tetF, tetT, tetM, H))
    result.push_back(cc.contours_used);

  return result;
}

int main(int argc, char *argv[]) {
  const char* trace_path = argv[1];
  const char* trace_name = argv[2];
  const int start = atoi(argv[3]);
  const int num_slices = atoi(argv[4]);

  SliceStack ss(trace_path, trace_name);

  // The first tile layer of contours.
  generated[start-1].push_back(new Tile(ss.getContoursAt(start)));

  for (int level = start; level < start + num_slices; level++) {
    if ((generated[level-1].size() == 0) || (ss.getSizeAt(level) <= 0))
      return -1;

    cout << "Level: " << level << endl;
    cout << "Contours: " << ss.getSizeAt(level) << endl;

    vector<set<int> > components = getHeatFlowValidComponents(ss, level);
    set<int> upper = ss.getContoursAt(level+1);

    // All nonunique tiles generated on this level.
    vector<Tile*> current_level_tiles;

    // Used to give each component the actual mesh.
    map<set<int>, vector<Component*> > contours_to_component;

    cout << "Generating tiles via exact set cover." << endl;

    // Go through all the previous level's valid tilings.
    for (Tile *parent : generated[level-1]) {
      bool last = (level == (start + num_slices - 1));

      // Do exact set cover and get all possible ways to cover.
      for (Tile *tile : Tiler::generateTiles(upper, parent, components, last)) {
        for (Component *component : tile->components)
          contours_to_component[component->contours_used].push_back(component);

        current_level_tiles.push_back(tile);
      }
    }

    cout << "Total tiles generated: " << current_level_tiles.size() << endl;

    // Generate mesh for each component.
    populateComponentMeshes(contours_to_component, ss);

    map<string, Tile*> best_tiles;

    cout << "Scoring each tile." << endl;

    // Pick the best of each connectivity.
    for (Tile *tile : current_level_tiles) {
      // Want to save the best from each unique parent connectivity.
      string tile_id = terribleHashFunction(tile);

      if (best_tiles.find(tile_id) == best_tiles.end()) {
        best_tiles[tile_id] = tile;
      }
      else {
        // Find the one with lower energy.
        double e1 = energy(tile);
        double e2 = energy(best_tiles[tile_id]);

        if (e1 < e2) {
          delete best_tiles[tile_id];
          best_tiles[tile_id] = tile;
        }
      }
    }

    // Populate the current level generated with the best tiles.
    for (auto map_it : best_tiles)
      generated[level].push_back(map_it.second);

    cout << "Valid tiles: " << generated[level].size() << endl;

    // Next bot is the current top.
    botV = topV;
    botF = topF;
    botM = topM;

    // For debugging.
    //   for (Tile *tile : generated[level])
    if (generated[level].size() > 0) {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      Eigen::VectorXi O;

      Tile *tile = generated[level].back();

      getTileMesh(tile, V, F, O, true);

      if ((level - start) % 5 == 0)
        viewTile(tile);

      igl::writeOFF("tile.off", V, F);
    }
  }

  cout << generated[start + num_slices - 1].size() << " valid tilings." << endl;

  // Show full tiles.
  for (Tile *tile : generated[start + num_slices - 1])
    viewTile(tile);
}
