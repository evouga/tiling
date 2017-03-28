#include <cstdio>

#include <igl/jet.h>
#include <igl/orient_outward.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "curvatureFlow.h"
#include "Helpers.h"

template<typename T, typename E>
int getOrig(const char* fn, const char* printf_str, E &O, T dummy) {
  FILE* inf = fopen(fn, "r");
  if (!inf) {
    perror("Could not open file");
    return -1;
  }
  printf("Getting lines...\n");

  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<T> orig;
  while((read = getline(&line, &len, inf)) != -1) {
    T val;
    sscanf(line, printf_str, &val);
    orig.push_back(val);
  }
  fclose(inf);

  printf("Found %lu elements\n", orig.size());
  O.resize(orig.size());
  for (int i = 0; i < orig.size(); ++i) {
    O(i) = orig[i];
  }

  return 0;
}

int main(int argc, char *argv[]) {
  if(argc < 2) {
    fprintf(stderr, "usage: %s input.[off,ply]\n",
            argv[0]);
    return -1;
  }

  Eigen::MatrixXd V, Vc;
  Eigen::MatrixXi F, Fc;
  Eigen::VectorXi C, Cc;
  bool orig_set = false;

  const char* inf = argv[1];
  if (strcmp(inf + strlen(inf) - 3, "off") == 0) {
    printf("Reading from input mesh of type OFF\n");
    igl::readOFF(inf, V,F);
    Eigen::VectorXi nums(F.rows()), other;
    nums.setZero();;
    igl::orient_outward(V,F,nums, F,other);

    // Look for something that looks like <input>_orig.txt for the original markers
    char* orig_fn = new char[strlen(inf) + 12];
    sprintf(orig_fn, "%.*s_orig.txt", (int)(strlen(inf) - 4), inf);
    printf("Looking for %s for original markers...\n", orig_fn);

    // Check to make sure this file exists.
    FILE* orig_f = fopen(orig_fn, "r");
    if (orig_f) {
      int i;
      getOrig(orig_fn, "%d", C, i);
      orig_set = true;
    } else {
      printf("Couldn't find original vertices; using dummy ones\n");
    }
  } else if (strcmp(inf + strlen(inf) - 3, "ply") == 0) {
    printf("Reading from input mesh of type PLY [%s]\n",
           inf + strlen(inf) - 3);
    igl::readPLY(inf, V,F);
  } else {
    printf("couldn't recogize file type of [%s] (type deduced to [%s]\n",
           inf, inf + strlen(inf) - 3);
    return -1;
  }

  if (!orig_set) {
    // Set the top, bottom, front, and back vertices to be static.
    Eigen::VectorXd mins = V.colwise().minCoeff();
    Eigen::VectorXd maxs = V.colwise().maxCoeff();

    C.resize(V.rows());
    C.setZero();
    /*
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
    */
    int fixed_idx = 0;
    // Set the middle line (x = bar x) to be constant
    Eigen::VectorXd avg = (mins + maxs) / 2.0;
    printf("avg is %lf", avg(fixed_idx));
    for (int i = 0; i < V.rows(); ++i) {
      if (std::abs(V(i, fixed_idx) - avg(fixed_idx)) < 1e-3) {
        C(i) = GLOBAL::original_marker;
      } else {
        C(i) = GLOBAL::nonoriginal_marker;
      }
    }

    bool do_extra = true;
    if (do_extra) {
      Eigen::VectorXi Cextra(C.rows());
      Cextra.setZero();
      // Let's make it two vertices
      for (int i = 0; i < F.rows(); ++i) {
        // If I've got a neighbor that's a marked vertex
        bool marked = false;
        Eigen::RowVector3d orig_r;
        for (int j = 0; j < F.cols(); ++j) {
          if (C(F(i, j)) == GLOBAL::original_marker) {
            marked = true;
            orig_r = V.row(F(i, j));
            break;
          }
        }
        // All vertices with X > marked are true.
        if (marked) {
          for (int j = 0; j < F.cols(); ++j) {
            int idx = F(i, j);
            if (V(idx, fixed_idx) > orig_r(fixed_idx)) {
              Cextra(idx) = GLOBAL::original_marker;
            }
          }
        }
      }
      // Add these to the mix
      for (int i = 0; i < C.rows(); ++i) {
        if (Cextra(i) != 0) {
          C(i) = GLOBAL::original_marker;
        }
      }
    }
  }
  igl::viewer::Viewer v;
  Eigen::MatrixXd cols;
  igl::jet(C, true, cols);
  v.data.clear();
  v.data.set_mesh(V, F);;
  v.data.set_colors(cols);
  v.launch();

  printf("Viewing before:\n");
  // Let's make sure it's manifold.
  Helpers::extractManifoldPatch(V, F, C);
  printf("Viewing after:\n");

  igl::jet(C, true, cols);
  v.data.clear();
  v.data.set_mesh(V, F);;
  v.data.set_colors(cols);
  v.launch();

  //biharmonic_view(V, F, C, Vc, true);
  double score = biharmonic_new(V, F, C, Vc, Fc, Cc);
  printf("Score: %lf\n", score);

  Eigen::MatrixXd tmp;
  biharmonic_view(V, F, C, tmp, false);

  //computeCurvatureFlow(V, F, C, 0.1, Vc);

  // igl::jet(Cc, true, cols);
  // v.data.clear();
  // v.data.set_mesh(Vc, Fc);
  // v.data.set_colors(cols);
  // v.launch();
}
