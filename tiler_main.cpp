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

int main(int argc, char *argv[]) {
  const char* trace_path = argv[1];
  const char* trace_name = argv[2];
  const int start = atoi(argv[3]);
  const int num_slices = atoi(argv[4]);

  SliceStack ss(trace_path, trace_name);

  for (int i = start; i < start + num_slices; i++) {
    int num_contours = ss.getSizeAt(i);
    cout << "Level: " << i << "\t# Contours: " << num_contours << endl;
    if (num_contours <= 0)
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

    TetMeshViewer::viewTetMesh(tetV, tetT, tetF, tetM, H);

    vector<ConnectedComponent> ccs = allPossibleTiles(tetV, tetF, tetT, tetM, H);

    for (const ConnectedComponent &cc : ccs) {
      Helpers::viewTriMesh(cc.V, cc.F, cc.M);
    }

    // The first tile layer of contours.
    // if (i == start) {
    //   generated[i].push_back(new Tile(upper));
    //   continue;
    // }

    // Go through all the previous level's valid tilings.
    // for (Tile *parent : generated[i-1]) {
    //   bool last = (i == (start + num_slices - 1));
    //   for (Tile *upper_tile : Tiler::generateTiles(upper, parent, last))
    //     generated[i].push_back(upper_tile);
    // }

    botV = topV;
    botF = topF;
    botM = topM;
  }

  // cout << endl << "Sample full tile stack top to bottom - " << endl;
  // cout << *generated[start + num_slices - 1][0] << endl;
}
