#include "Helpers.h"

#include "curvatureFlow.h"
#include "glob_defs.h"

#include <iostream>
#include <cstdio>
#include <set>

#include <Eigen/Core>

#include <igl/collapse_small_triangles.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/cotmatrix.h>
#include <igl/jet.h>
#include <igl/remove_duplicates.h>
#include <igl/remove_unreferenced.h>
#include <igl/resolve_duplicated_faces.h>
#include <igl/triangle/triangulate.h>
#include <igl/viewer/Viewer.h>
#include <igl/writeOFF.h>

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
  sprintf(buffer, "Qpq%fY", GLOBAL::TET_RATIO);

  // This prints out a number for some reason.
  printf("[%s:%u] The igl tetrahedralize function produces this output: ", __FILE__, __LINE__);
  igl::copyleft::tetgen::tetrahedralize(V, F, VM, FM, buffer, TV, TT, TF, TO);
  printf("\n");

  // Tetgen has some weird markers it introduces. Make them all either
  // original or nonoriginal.
  for (int i = 0; i < TO.rows(); ++i)
    TO(i) = max(TO(i), GLOBAL::nonoriginal_marker);
}

void triangulate(const Eigen::MatrixXd &P, const Eigen::MatrixXi &E,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  char buffer[80];
  sprintf(buffer, "QDYa%fq", GLOBAL::TRI_AREA);

  // For holes, not used right now.
  Eigen::MatrixXd H_unused(0, 2);

  igl::triangle::triangulate(P, E, H_unused, buffer, V, F);
}


namespace {

void resetView (igl::viewer::Viewer &view,
                const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &M, const Eigen::MatrixXd &C) {
  view.data.clear();

  view.data.set_mesh(V, F);
  view.data.set_colors(C);

  // Show only one marker for each contour.
  set<int> used;

  for (int i = 0; i < M.rows(); i++) {
    if (M(i) != GLOBAL::nonoriginal_marker && used.find(M(i)) == used.end()) {
      view.data.add_label(V.row(i), to_string(M(i)));
      used.insert(M(i));
    }
  }
}

}

void viewTriMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 const Eigen::VectorXi &M) {
  igl::viewer::Viewer v;

  // Colors.
  Eigen::MatrixXd C;
  igl::jet(M, true, C);

  Eigen::MatrixXd V_biharmonic = V;
  Eigen::MatrixXi F_biharmonic = F;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);

  v.callback_key_down = [&](igl::viewer::Viewer& v, unsigned char key, int modifier) {
    if (key == ' ') {
      biharmonic(V_biharmonic, F_biharmonic, M, L, V_biharmonic);
      //igl::collapse_small_triangles(V_biharmonic, F, 1e-6, F_biharmonic);
      printf("Finished biharmonic\n");
      resetView(v, V_biharmonic, F_biharmonic, M, C);
    }
    else if (key == 'R') {
      printf("Resetting matrices\n");
      igl::cotmatrix(V, F, L);
      V_biharmonic = V;
      F_biharmonic = F;
      resetView(v, V_biharmonic, F_biharmonic, M, C);
    } else if (key == 'B') {
      printf("Resetting laplacian matrix\n");
      igl::cotmatrix(V,F,L);
    }
    return true;
  };

  resetView(v, V, F, M, C);

  // Getting random segfaults. Write it to see if it is a mesh problem.
  igl::writeOFF("mesh.off", V, F);

  printf("  Press ' ' (space) to compute biharmonic\n");
  printf("  Press 'R' to reset the vertices back to original\n");
  printf("  Press 'B' to re-compute Laplacian matrix\n");

  v.launch();
}

void combineMesh(const vector<Eigen::MatrixXd> &Vs,
                 const vector<Eigen::MatrixXi> &Fs,
                 const vector<Eigen::VectorXi> &Ms,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &M) {
  assert(Vs.size() == Fs.size());
  assert(Vs.size() == Ms.size());

  int total_vertices = 0;
  int total_faces = 0;

  for (const Eigen::MatrixXd &vertices : Vs)
    total_vertices += vertices.rows();
  for (const Eigen::MatrixXi &faces : Fs)
    total_faces += faces.rows();

  V.resize(total_vertices, 3);
  F.resize(total_faces, 3);
  M.resize(total_vertices);

  int vertex_offset = 0;
  int face_offset = 0;

  for (int i = 0; i < Vs.size(); i++) {
    const Eigen::MatrixXd &vertices = Vs[i];
    const Eigen::MatrixXi &faces = Fs[i];
    const Eigen::VectorXi &markers = Ms[i];

    for (int j = 0; j < vertices.rows(); j++)
      V.row(vertex_offset + j) = vertices.row(j);

    for (int j = 0; j < faces.rows(); j++) {
      for (int k = 0; k < 3; k++)
        F(face_offset + j, k) = vertex_offset + faces(j, k);
    }

    for (int j = 0; j < markers.rows(); j++)
      M(vertex_offset + j) = markers(j);

    vertex_offset += vertices.rows();
    face_offset += faces.rows();
  }

  removeDuplicates(V, F, M);
}

void extractShell(Eigen::MatrixXd &V1, Eigen::MatrixXi &F1, Eigen::VectorXi &M1,
                  Eigen::MatrixXd &V2, Eigen::MatrixXi &F2, Eigen::VectorXi &M2) {
  Eigen::MatrixXi F_unique;
  Eigen::VectorXi J;
  igl::resolve_duplicated_faces(F1, F_unique, J);

  // J is the same size as V
  igl::remove_unreferenced(V1, F_unique, V2, F2, J);

  M2.resize(V2.rows());

  // J_i is -1 if the vertex was removed.
  for (int i = 0; i < M1.rows(); ++i) {
    if (J(i) != -1)
      M2(J(i)) = M1(i);
  }
}

} // namespace Helpers
