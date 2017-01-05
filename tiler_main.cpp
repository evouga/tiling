#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <cstdlib>

#include <Eigen/Core>

#include "Tiler.h"
#include "SliceStack.h"
#include "glob_defs.h"
#include "viewTetMesh.h"
#include "TilingUtils.h"
#include "Helpers.h"

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

int main(int argc, char *argv[]) {
  const char* trace_path = argv[1];
  const char* trace_name = argv[2];
  const int start = atoi(argv[3]);
  const int num_slices = atoi(argv[4]);

  SliceStack ss(trace_path, trace_name);

  // The first tile layer of contours.
  generated[start-1].push_back(new Tile(ss.getContoursAt(start)));

  for (int i = start; i < start + num_slices; i++) {
    if (ss.getSizeAt(i) <= 0)
      return -1;

    // Need heat values for the contour tree.
    ss.triangulateSlice(i, GLOBAL::TRI_AREA,
                        botV, botF, topV, topF, botM, topM);

    // Helpers::viewTriMesh(botV, botF, botM);
    // Helpers::viewTriMesh(topV, topF, topM);

    ss.tetrahedralizeSlice(botV, botF, topV, topF,
                           botM, topM, tetV, tetT, tetF, tetM);

    // TetMeshViewer::viewTetMesh(tetV, tetT, tetF, tetM, H);

    ss.computeLaplace(tetV, tetT, tetF, tetM, H);

    // TetMeshViewer::viewTetMesh(tetV, tetT, tetF, tetM, H);

    vector<ConnectedComponent> ccs = allPossibleTiles(tetV, tetF, tetT, tetM, H);
    vector<set<int> > components = convertCCsToCover(ccs);

    // for (const ConnectedComponent &cc : ccs) {
    //   Helpers::viewTriMesh(cc.V, cc.F, cc.M);
    // }

    set<int> upper = ss.getContoursAt(i+1);

    // Go through all the previous level's valid tilings.
    for (Tile *parent : generated[i-1]) {
      bool last = (i == (start + num_slices - 1));
      for (Tile *tile : Tiler::generateTiles(upper, parent, components, last)) {

        // TODO(bradyz): Combine both component classes.
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

        generated[i].push_back(tile);
      }
    }

    cout << "Level: " << i << endl;
    cout << "Contours: " << ss.getSizeAt(i) << endl;
    cout << "Valid tiles: " << generated[i].size() << " tilings." << endl;

    // Next bot is the current top.
    botV = topV;
    botF = topF;
    botM = topM;
  }

  cout << generated[start + num_slices - 1].size() << " valid tilings." << endl;

  // Show full tiles.
  for (Tile *tile : generated[start + num_slices - 1]) {
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
    Helpers::viewTriMesh(tileV, tileF, tileM);
  }
}
