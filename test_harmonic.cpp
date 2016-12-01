#include <cstdio>

#include <igl/readOFF.h>
#include <igl/readPLY.h>

#include "glob_defs.h"
#include "curvatureFlow.h"

int main(int argc, char *argv[]) {
  if(argc < 2) {
    fprintf(stderr, "usage: %s input.[off,ply]\n",
            argv[0]);
    return -1;
  }

  Eigen::MatrixXd V, Vc;
  Eigen::MatrixXi F;
  Eigen::VectorXi C;

  const char* inf = argv[1];
  if (strcmp(inf + strlen(inf) - 3, "off") == 0) {
    printf("Reading from input mesh of type OFF\n");
    igl::readOFF(inf, V,F);
  } else if (strcmp(inf + strlen(inf) - 3, "ply") == 0) {
    printf("Reading from input mesh of type PLY [%s]\n",
           inf + strlen(inf) - 3);
    igl::readPLY(inf, V,F);
  } else {
    printf("couldn't recogize file type of [%s] (type deduced to [%s]\n",
           inf, inf + strlen(inf) - 3);
    return -1;
  }

  // Set the top, bottom, front, and back vertices to be static.
  Eigen::VectorXd mins = V.colwise().minCoeff();
  Eigen::VectorXd maxs = V.colwise().maxCoeff();

  C.resize(V.rows());
  C.setZero();
  for (int i = 0; i < V.rows(); ++i) {
    bool original = false;
    for (int j = 0; j < 3; ++j) {
      if (V(i, j) < mins(j)*0.95 || V(i, j) > maxs(j)*0.95) {
        original = true;
      }
    }

    if (original) {
      C(i) = GLOBAL::original_marker;
    } else {
      C(i) = GLOBAL::nonoriginal_marker;
    }
  }

  biharmonic(V, F, C, Vc);
  //computeCurvatureFlow(V, F, C, 0.1, Vc);
}
