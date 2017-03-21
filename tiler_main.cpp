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
string terribleHashFunction (Tile *tile) {
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

void combineComponentsIntoMesh(const vector<Component*> &components,
                               Eigen::MatrixXd &V,
                               Eigen::MatrixXi &F,
                               Eigen::VectorXi &O) {
  vector<Eigen::MatrixXd> tileVs;
  vector<Eigen::MatrixXi> tileFs;
  vector<Eigen::VectorXi> tileMs;

  for (Component *component : components) {
    tileVs.push_back(component->V);
    tileFs.push_back(component->F);
    tileMs.push_back(component->M);
  }

  Eigen::MatrixXd tileV;
  Eigen::MatrixXi tileF;
  Eigen::VectorXi tileM;

  Helpers::combineMesh(tileVs, tileFs, tileMs, tileV, tileF, tileM);
  Helpers::extractShell(tileV, tileF, tileM, V, F, O);
}

double energy(Tile *tile) {
  double score = 0.0;

  vector<Eigen::MatrixXd> tileVs;
  vector<Eigen::MatrixXi> tileFs;
  vector<Eigen::VectorXi> tileMs;

  vector<Component*> components = getTileComponents(tile, 1);

  for (const Component *component : components) {
    const Eigen::MatrixXd &V = component->V;
    const Eigen::MatrixXi &F = component->F;
    const Eigen::VectorXi &O = component->M;

    Eigen::MatrixXd V_shell, V_biharmonic;
    Eigen::MatrixXi F_shell, F_biharmonic;
    Eigen::VectorXi O_shell, O_biharmonic;

    // Make sure to pass in a triangle mesh.
    Helpers::extractShell(V, F, O, V_shell, F_shell, O_shell);
    score += biharmonic_new(V_shell, F_shell, O_shell,
                            V_biharmonic, F_biharmonic, O_biharmonic);
  }

  Eigen::MatrixXd V, V_unused;
  Eigen::MatrixXi F, F_unused;
  Eigen::VectorXi M, M_unused;

  // Why are these not equal? :(.
  combineComponentsIntoMesh(getTileComponents(tile, 1), V, F, M);
  cout << "energy: " << score << endl;
  cout << "together: " << biharmonic_new(V, F, M,
                                         V_unused, F_unused, M_unused) << endl;

  // TODO: add genus computation per connected component.
  // double genus = (V.rows() - 0.5 * F.rows() - 2.0) / -2.0;

  if (std::isnan(score))
    return INFINITY;

  return score;
}

void viewTile(Tile *tile, int num_tiles=-1, bool save=false) {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi M;

  combineComponentsIntoMesh(getTileComponents(tile, num_tiles), V, F, M);
  Helpers::viewTriMesh(V, F, M);

  if (save) {
    igl::writeOFF("tile.off", V, F);

    // Also write the original vertices
    FILE* of = fopen("tile_orig.txt", "w");
    for (int i = 0; i < M.rows(); ++i)
      fprintf(of, "%d\n", M(i));
    fclose(of);
  }
}

map<set<int>, ConnectedComponent> getHeatFlowValidComponents(SliceStack &ss,
                                                             int level) {
  ss.triangulateSlice(level, GLOBAL::TRI_AREA,
                      botV, botF, topV, topF, botM, topM);
  ss.tetrahedralizeSlice(botV, botF, topV, topF,
                         botM, topM, tetV, tetT, tetF, tetM);
  ss.computeLaplace(tetV, tetT, tetF, tetM, H);

  return possibleTileMap(tetV, tetF, tetT, tetM, H);
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

    // Map from contours used to ConnectedComponent mesh.
    map<set<int>, ConnectedComponent> contours_to_component =
        getHeatFlowValidComponents(ss, level);

    vector<set<int> > components;
    for (auto it : contours_to_component)
      components.push_back(it.first);

    set<int> upper = ss.getContoursAt(level+1);

    // All nonunique tiles generated on this level.
    vector<Tile*> current_level_tiles;

    cout << "Generating tiles via exact set cover." << endl;

    // Go through all the previous level's valid tilings.
    for (Tile *parent : generated[level-1]) {
      bool last = (level == (start + num_slices - 1));

      // Do exact set cover and get all possible ways to cover.
      for (Tile *tile : Tiler::generateTiles(upper, parent, components, last)) {

        // Go through component and attach the mesh.
        for (Component *component : tile->components) {
          component->V = contours_to_component[component->contours_used].V;
          component->F = contours_to_component[component->contours_used].F;
          component->M = contours_to_component[component->contours_used].M;
        }

        current_level_tiles.push_back(tile);
      }
    }

    cout << "Total tiles generated: " << current_level_tiles.size() << endl;

    map<string, Tile*> best_tiles;

    cout << "Scoring each tile." << endl;

    // Pick the best of each connectivity.
    for (Tile *tile : current_level_tiles) {
      // Want to save the best from each unique parent connectivity.
      string tile_id = terribleHashFunction(tile);
      printf("Tile:\n=====\n%s=====\n", tile_id.c_str());
      printf("Energy is %lf\n", energy(tile));

      if (best_tiles.find(tile_id) == best_tiles.end()) {
        printf("Empty, adding to slot\n");
        best_tiles[tile_id] = tile;
      }
      else {
        // Find the one with lower energy.
        double e1 = energy(tile);
        double e2 = energy(best_tiles[tile_id]);

        printf("First Score: %lf\n", e1);
        viewTile(tile, 1);
        printf("Second Score: %lf\n", e2);
        viewTile(best_tiles[tile_id], 1);

        if (e1 < e2) {
          // Purge the loser.
          delete best_tiles[tile_id];
          best_tiles[tile_id] = tile;
        }
      }
    }

    // Populate the current level generated with the best tiles.
    for (auto map_it : best_tiles) {
      generated[level].push_back(map_it.second);
    }

    cout << "Valid tiles: " << generated[level].size() << endl;

    // Next bot is the current top.
    botV = topV;
    botF = topF;
    botM = topM;

    // For debugging.
    if (generated[level].size() > 0) {
      if (level > start &&
          ((level - start) % 5 == 0 || (level - start == num_slices - 1))) {
        for (Tile *tile : generated[level])
          viewTile(tile);

        cout << "Saving a full tile for debugging." << endl;
        viewTile(generated[level].back(), -1, true);
      }
    }
  }

  cout << generated[start + num_slices - 1].size() << " valid tilings." << endl;

  // Show full tiles.
  for (Tile *tile : generated[start + num_slices - 1])
    viewTile(tile);
}
