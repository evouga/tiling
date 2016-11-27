#include "Helpers.h"
#include "glob_defs.h"

#include <iostream>
#include <cstdio>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/remove_duplicates.h>
#include <igl/triangle/triangulate.h>

using namespace std;

namespace Helpers {

void removeDuplicates(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                      Eigen::VectorXi &M) {
  // Make sure the markers correspond to the vertices.
  assert(V.rows() == M.rows());

  Eigen::MatrixXd fixedV;
  Eigen::MatrixXi fixedF;
  Eigen::VectorXi toUnique;
  igl::remove_duplicates(V, F, fixedV, fixedF, toUnique, GLOBAL::EPS);

  // Populate the new markers.
  Eigen::VectorXi fixedM(fixedV.rows());
  fixedM.setZero();
  for (int old_i = 0; old_i < M.rows(); ++old_i) {
    int new_i = toUnique(old_i);
    // A non-original marker can be merged into an original one.
    fixedM(new_i) = max(fixedM(new_i), M(old_i));
  }

  if (GLOBAL::DEBUG) {
    cout << "Removed duplicate vertices." << endl;
    cout << "V: " << V.rows() << " -> " << fixedV.rows() << endl;
    cout << "F: " << F.rows() << " -> " << fixedF.rows() << endl;
  }

  // Make the method seem in-place.
  V = fixedV;
  F = fixedF;
  M = fixedM;
}

Eigen::VectorXi getFaceMarkers(const Eigen::MatrixXi &F,
                               const Eigen::VectorXi &M) {
  Eigen::VectorXi FM(F.rows());
  for (int i = 0; i < F.rows(); ++i) {
    bool original = (M(F(i, 0)) == GLOBAL::original_marker &&
                     M(F(i, 1)) == GLOBAL::original_marker &&
                     M(F(i, 2)) == GLOBAL::original_marker);
    FM(i) = (original) ? (GLOBAL::original_marker) : (GLOBAL::nonoriginal_marker);
  }
  return FM;
}

void tetrahedralize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    const Eigen::VectorXi &VM, const Eigen::VectorXi &FM,
                    Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                    Eigen::MatrixXi &TF, Eigen::VectorXi &TO) {
  // Make sure markers are populated for vertices and faces.
  assert(V.rows() == VM.rows() && F.rows() == FM.rows());

  char buffer[80];
  sprintf(buffer, "pq%fY", GLOBAL::TET_RATIO);

  igl::copyleft::tetgen::tetrahedralize(V, F, VM, FM, buffer, TV, TT, TF, TO);

  if (GLOBAL::DEBUG) {
    cout << "Tetrahedralized mesh." << endl;
    cout << "F: " << F.rows() << " -> " << TF.rows() << endl;
  }
}

void triangulate(const Eigen::MatrixXd &P, const Eigen::MatrixXi &E,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  char buffer[80];
  sprintf(buffer, "DYa%fq", GLOBAL::TRI_AREA);

  // For holes, not used right now.
  Eigen::MatrixXd H_unused(0, 2);

  igl::triangle::triangulate(P, E, H_unused, buffer, V, F);
}

} // namespace Helpers
