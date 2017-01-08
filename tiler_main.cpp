#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <Eigen/Core>

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

vector<set<int> > convertCCsToCover(const vector<ConnectedComponent> &ccs) {
  vector<set<int> > result;
  for (const ConnectedComponent &cc : ccs)
    result.push_back(cc.contours_used);
  return result;
}

string terribleHashFunction (const Tile *tile) {
  vector<set<int> > connected = tile->getUpperConnected();
  sort(connected.begin(), connected.end());

  stringstream ss;
  for (set<int> foo : connected) {
    for (int bar : foo)
      ss << bar << " ";
    ss << endl;
  }
  return ss.str();
}

// Extracts the surface of the current level of the given tile.
void getTileMesh (Tile *tile,
                  Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O) {
  vector<Eigen::MatrixXd> tileVs;
  vector<Eigen::MatrixXi> tileFs;
  vector<Eigen::VectorXi> tileMs;

  for (Component *component : tile->components) {
    tileVs.push_back(component->V);
    tileFs.push_back(component->F);
    tileMs.push_back(component->M);
  }

  // Full mesh.
  Eigen::MatrixXd tileV;
  Eigen::MatrixXi tileF;
  Eigen::VectorXi tileM;

  Helpers::combineMesh(tileVs, tileFs, tileMs, tileV, tileF, tileM);
  Helpers::extractShell(tileV, tileF, tileM, V, F, O);
}

double energy (Tile *tile) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi O;
  getTileMesh(tile, V, F, O);

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
  cout << endl << "Sample full tile stack top to bottom - " << endl;
  cout << *tile << endl;

  // Get vector of full mesh.
  vector<Component*> components = tile->getAllComponents();

  vector<Eigen::MatrixXd> tileVs;
  vector<Eigen::MatrixXi> tileFs;
  vector<Eigen::VectorXi> tileMs;

  for (Component *component : components) {
    tileVs.push_back(component->V);
    tileFs.push_back(component->F);
    tileMs.push_back(component->M);
  }

  // Full mesh.
  Eigen::MatrixXd tileV;
  Eigen::MatrixXi tileF;
  Eigen::VectorXi tileM;

  Helpers::combineMesh(tileVs, tileFs, tileMs, tileV, tileF, tileM);

  Eigen::MatrixXd Vborder;
  Eigen::MatrixXi Fborder;
  Eigen::VectorXi Oborder;

  Helpers::extractShell(tileV, tileF, tileM, Vborder, Fborder, Oborder);
  Helpers::viewTriMesh(Vborder, Fborder, Oborder);
}

int main(int argc, char *argv[]) {
  const char* trace_path = argv[1];
  const char* trace_name = argv[2];
  const int start = atoi(argv[3]);
  const int num_slices = atoi(argv[4]);

  SliceStack ss(trace_path, trace_name);

  // The first tile layer of contours.
  generated[start-1].push_back(new Tile(ss.getContoursAt(start)));

  for (int i = start; i < start + num_slices; i++) {
    if ((generated[i-1].size() == 0) || (ss.getSizeAt(i) <= 0))
      return -1;

    // Need heat values for the contour tree.
    ss.triangulateSlice(i, GLOBAL::TRI_AREA,
                        botV, botF, topV, topF, botM, topM);
    ss.tetrahedralizeSlice(botV, botF, topV, topF,
                           botM, topM, tetV, tetT, tetF, tetM);
    ss.computeLaplace(tetV, tetT, tetF, tetM, H);

    vector<ConnectedComponent> ccs = allPossibleTiles(tetV, tetF, tetT, tetM, H);
    vector<set<int> > components = convertCCsToCover(ccs);

    set<int> upper = ss.getContoursAt(i+1);

    vector<Tile*> current_level_tiles;

    // Go through all the previous level's valid tilings.
    for (Tile *parent : generated[i-1]) {
      bool last = (i == (start + num_slices - 1));
      for (Tile *tile : Tiler::generateTiles(upper, parent, components, last))
        current_level_tiles.push_back(tile);
    }

    map<string, Tile*> best_tiles;

    // Pick the best of each connectivity.
    for (Tile *tile : current_level_tiles) {
      // TODO(bradyz): Combine both component classes.
      // NOTE: This loop must be done before calling getTileMesh, energy.
      for (Component *component : tile->components) {
        int index = -1;
        for (int j = 0; j < components.size(); j++) {
          if (component->contours_used == components[j])
            index = j;
        }
        component->V = ccs[index].V;
        component->F = ccs[index].F;
        component->M = ccs[index].M;
      }

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
      generated[i].push_back(map_it.second);

    cout << "Level: " << i << endl;
    cout << "Contours: " << ss.getSizeAt(i) << endl;
    cout << "Valid tiles: " << generated[i].size() << endl;

    // Next bot is the current top.
    botV = topV;
    botF = topF;
    botM = topM;
  }

  cout << generated[start + num_slices - 1].size() << " valid tilings." << endl;

  // Show full tiles.
  for (Tile *tile : generated[start + num_slices - 1]) {
    viewTile(tile);
  }
}
