#define MESHFIX_WITH_EIGEN

#include "Helpers.h"

#include "curvatureFlow.h"
#include "glob_defs.h"
#include "marching_tets.h"

//#include "meshfix.h"

#include <iostream>
#include <cstdio>
#include <set>
#include <utility>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/cotmatrix.h>
#include <igl/extract_manifold_patches.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/jet.h>
#include <igl/remove_duplicates.h>
#include <igl/remove_unreferenced.h>
#include <igl/resolve_duplicated_faces.h>
#include <igl/triangle/triangulate.h>
#include <igl/viewer/Viewer.h>
#include <igl/writeOFF.h>
#include <igl/collapse_small_triangles.h>

using namespace std;

namespace Helpers {

// Determines if each edge is manifold or not.
bool is_edge_manifold(const Eigen::MatrixXi &F,
                      Eigen::VectorXi &M) {
  map<pair<int,int>, int> edge_count;
  vector<vector<pair<int,int> > > f2e;

  // Get all the edge counts, and the map from edges back to faces.
  for (int i = 0; i < F.rows(); ++i) {
    f2e.push_back({});
    for (int j = 0; j < F.cols(); ++j) {
      int next_j = (j + 1) % 3;

      int u = F(i, j);
      int v = F(i, next_j);

      if (u > v) {int t = u; u = v; v = t;}

      edge_count[make_pair(u, v)] ++;
      f2e[i].push_back(make_pair(u, v));
    }
  }

  bool is_manifold = true;

  // Report which ones are non-manifold.
  M.resize(F.rows());
  M.setZero();
  for (int i = 0; i < F.rows(); ++i) {
    for (const auto& e : f2e[i]) {
      if (edge_count[e] > 2) {
        M(i) = max(M(i), edge_count[e]);
        is_manifold = false;
      }
    }
  }

  return is_manifold;
}


bool isManifold(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXi &M,
    bool changeMarkers) {
  bool is_manifold = igl::is_edge_manifold(F);
  Eigen::VectorXi B;
  bool is_v_manifold = igl::is_vertex_manifold(F, B);
  if (!is_v_manifold) {
    printf("Warning: Your mesh is not vertex manifold!\n");
    for (int i = 0; i < M.rows(); ++i) {
      if (!B(i)) {
        printf("Vertex %d has %d (was %d)\n", i, B(i), M(i));
        if (changeMarkers)
          M(i) = 100;
      }
    }
  }
  if (!is_manifold) {
    printf("Warning: Your mesh is not edge manifold!!\n");
    Eigen::VectorXi NME;
    is_edge_manifold(F, NME);
    printf("Errors at: \n");
    for (int i = 0; i < F.rows(); ++i) {
      if (NME(i)) {
        printf("  row %d repeated %d\n", i, NME(i));
        for (int j = 0; j < F.cols(); ++j) {
          if (changeMarkers)
            M(F(i,j)) = 100;
          printf("   -> Setting %d to 100\n", F(i,j));
        }
      }
    }
  }
  return is_manifold && is_v_manifold;
}
// Overwrites inputs.
// Each manifold patch must have at least minFaces faces or it will be
// removed.
bool extractManifoldPatch(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                     Eigen::VectorXi &M, int minFaces, bool changeMarkers) {

  if(isManifold(V, F, M, changeMarkers)) return true;

  Eigen::VectorXi patchId;
  int nPatch = igl::extract_manifold_patches(F, patchId);
  // Find the largest one.
  int* patchCount = new int[nPatch];
  // Set to be zero.
  for (int i = 0; i < nPatch; ++i) {
    patchCount[i] = 0;
  }
  // Find the counts of patches.
  for (int i = 0; i < patchId.rows(); ++i) {
    patchCount[patchId(i)]++;
  }
  // Get the total number of faces.
  int totalRows = 0;
  for (int i = 0; i < nPatch; ++i) {
    if (patchCount[i] >= minFaces) {
      totalRows += patchCount[i];
    }
  }

  printf("Removing %d faces!\n", F.rows() - totalRows);
  // Remove all faces that are not associated with maxIdx.
  Eigen::MatrixXi newF(totalRows, F.cols());
  int newIdx = 0;
  for (int i = 0; i < F.rows(); ++i) {
    // Add it to the new mesh if it belongs.
    int pid = patchId(i);
    if (patchCount[pid] >= minFaces) {
      newF.row(newIdx++) = F.row(i);
    }
  }

  // Then, remove unreferenced vertices and update markers.
  Eigen::MatrixXd V2;
  Eigen::VectorXi M2, J;
  // J is indices into V
  // Just write directly into F.
  igl::remove_unreferenced(V, newF, V2, F, J);

  M2.resize(V2.rows());

  // J_i is -1 if the vertex was removed.
  for (int i = 0; i < M.rows(); ++i) {
    if (J(i) != -1)
      M2(J(i)) = M(i);
  }

  V = V2;
  M = M2;

  delete[] patchCount;

  return false;
}

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
    if (fixedM(new_i) == 0 || fixedM(new_i) == GLOBAL::nonoriginal_marker) {
      fixedM(new_i) = M(old_i);
    } else {
      // If the previous is not a non-original marker, then just set it to
      // be the minimum of the previous two.
      fixedM(new_i) = min(fixedM(new_i), M(old_i));
    }
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
  printf("[%s:%u] igl:tetrahedralize: ", __FILE__, __LINE__);
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

void set_viewer(igl::viewer::Viewer &viewer,
                const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &M) {
  if (!isMeshOkay(V, F))
    cout << "Mesh has problems." << endl;

  // Add colors.
  Eigen::MatrixXd C;
  igl::jet(M, true, C);

  viewer.data.clear();
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(C);

  // Show only one marker for each contour.
  set<int> used;

  for (int i = 0; i < M.rows(); i++) {
    if (M(i) == GLOBAL::nonoriginal_marker)
      continue;
    else if (used.find(M(i)) != used.end())
      continue;

    viewer.data.add_label(V.row(i), to_string(M(i)));

    if (M(i) == -100 || M(i) == 100)
      continue;

    used.insert(M(i));
  }
}

void set_viewer_with_color(igl::viewer::Viewer &viewer,
                           const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::VectorXd &C_unnormalized) {
  // Take the log of the squared values.
  Eigen::VectorXd C_normalized(C_unnormalized.rows());
  for (int i = 0; i < C_unnormalized.rows(); i++)
    C_normalized(i) = log(C_unnormalized(i) * C_unnormalized(i) + 1e-5);

  // Turn scalars into colors.
  Eigen::MatrixXd C;
  igl::jet(C_normalized, true, C);

  viewer.data.clear();
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(C);
}

}

void viewTriMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 const Eigen::VectorXi &O) {
  igl::viewer::Viewer viewer;

  Eigen::MatrixXd V_biharmonic = V;
  Eigen::MatrixXi F_biharmonic = F;
  Eigen::VectorXi O_biharmonic = O;
  Eigen::SparseMatrix<double> L;
  vector<int> new_vertices;

  int times_called = 0;

  viewer.callback_key_down = [&](igl::viewer::Viewer& viewer,
                                 unsigned char key, int modifier) {
    if (key == ' ') {
      Eigen::MatrixXd V_next = V_biharmonic;
      Eigen::MatrixXi F_next = F_biharmonic;
      Eigen::VectorXi M_next = O_biharmonic;

      double energy = 0.0;

      if (times_called == 0) {
        energy = biharmonic_new(V_biharmonic, F_biharmonic, O_biharmonic,
                                V_next, F_next, M_next, &new_vertices);
        // Compute L on the mesh now.
        igl::cotmatrix(V_next, F_next, L);
      }
      else {
        energy = biharmonic(V_biharmonic, F_biharmonic,
                            O_biharmonic, O_biharmonic, // same fixed verts.
                            L,
                            V_next);
      }

      V_biharmonic = V_next;
      F_biharmonic = F_next;
      O_biharmonic = M_next;

      printf("Finished biharmonic: %lf\n", energy);
      set_viewer(viewer, V_biharmonic, F_biharmonic, O_biharmonic);

      times_called++;
    }
    else if (key == 'R') {
      printf("Resetting mesh.\n");
      V_biharmonic = V;
      F_biharmonic = F;
      O_biharmonic = O;

      times_called = 0;
      set_viewer(viewer, V_biharmonic, F_biharmonic, O_biharmonic);
    } else if (key == 'L') {
      printf("Recomputing Laplacian\n");
      // Compute L on the mesh now.
      igl::cotmatrix(V_biharmonic, F_biharmonic, L);
    } else if (key == 'H') {
      printf("Viewing errors...\n");
      std::set<int> vissues;
      isMeshOkay(V_biharmonic, F_biharmonic, vissues);
      // Set the vertices with issues.
      for (int i : vissues) {
        O_biharmonic(i) = 100;
      }
      set_viewer(viewer, V_biharmonic, F_biharmonic, O_biharmonic);
    } else if (key == 'C') {
      printf("Visualizing energy at each vertex.\n");

      // All the new points should be ignored in energy calculation.
      Eigen::VectorXd energy_density = biharmonic_energy_per_vertex(V_biharmonic,
                                                                    F_biharmonic,
                                                                    new_vertices);

      set_viewer_with_color(viewer, V_biharmonic, F_biharmonic, energy_density);
    /*} else if (key == 'F') {
      Eigen::MatrixXd W;
      Eigen::MatrixXi G;

      //meshfix(V_biharmonic, F_biharmonic, W, G);

      Eigen::VectorXi O_unused(W.rows());
      set_viewer(viewer, W, G, O_unused);
      */
    }

    return true;
  };

  set_viewer(viewer, V, F, O);

  // Getting random segfaults. Write it to see if it is a mesh problem.
  igl::writeOFF("mesh.off", V, F);

  printf("  Press ' ' (space) to compute biharmonic\n");
  printf("  Press 'R' to reset the vertices back to original\n");

  viewer.launch();
}

void combineMesh(const vector<Eigen::MatrixXd> &Vs,
                 const vector<Eigen::MatrixXi> &Fs,
                 const vector<Eigen::VectorXi> &Ms,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &M) {
  // Vectors should all be the same length.
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

void extractShell(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1,
                  const Eigen::VectorXi &M1,
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

  // Also make sure it's manifold.
  extractManifoldPatch(V2, F2, M2);
}

void collapseSmallTriangles(const Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                            double eps) {
  Eigen::MatrixXi F_tmp;
  igl::collapse_small_triangles(V, F, eps, F_tmp);

  cout << F.rows() - F_tmp.rows() << " faces removed." << endl;
  F = F_tmp;
}

// Remove unused vertices.
void removeUnreferenced(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                        Eigen::VectorXi &O) {
  int original_size = V.rows();

  Eigen::MatrixXd V_new;
  Eigen::MatrixXi F_new;
  Eigen::VectorXi old_to_new;

  igl::remove_unreferenced(V, F, V_new, F_new, old_to_new);

  Eigen::VectorXi O_new(V_new.rows());

  for (int old_index = 0; old_index < old_to_new.rows(); old_index++) {
    int new_index = old_to_new(old_index);

    // New index is -1 if it was removed..
    if (new_index != -1)
      O_new(new_index) = O(old_index);
  }

  if (GLOBAL::DEBUG)
    cout << (V.rows() - V_new.rows()) << " unused vertices removed." << endl;

  // Make it in place.
  V = V_new;
  F = F_new;
  O = O_new;
}

bool isMeshOkay(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double eps) {
  std::set<int> ignore;
  return isMeshOkay(V, F, ignore, eps);
}

bool isMeshOkay(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                std::set<int> &vertex_issues, double eps) {
  bool result = true;
  double min_area = 1e10;

  set<int> vertices;

  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      // Check if vertex indices are valid.
      if (F(i, j) < 0 || F(i, j) >= V.rows()) {
        cout << "Vertex index " << F(i, j) << " is out of bounds." << endl;
        vertex_issues.insert(F(i, j));
        result = false;
      }

      vertices.insert(F(i, j));
    }

    // Check if faces have small areas.
    Eigen::RowVector3d u = V.row(F(i, 1)) - V.row(F(i, 0));
    Eigen::RowVector3d v = V.row(F(i, 2)) - V.row(F(i, 0));
    double area = u.cross(v).norm();

    if (area < eps) {
      cout << "Face " << i << " has area " << area << "("
           << F(i, 0) << ", " << F(i, 1) << ", " << F(i, 2) << ")\n";
      result = false;
      vertex_issues.insert(F(i, 0));
      vertex_issues.insert(F(i, 1));
      vertex_issues.insert(F(i, 2));
    }

    min_area = min(min_area, area);
  }

  int unused = 0;

  // Check if all vertex indices are used.
  for (int i = 0; i < V.rows(); i++) {
    if (vertices.find(i) == vertices.end()) {
      unused++;

      cout << "Vertex " << i << " not used in faces." << endl;
      result = false;
      vertex_issues.insert(i);
    }
  }

  // Check for nans in vertex positions.
  for (int i = 0; i < V.rows(); i++) {
    bool has_nan = false;

    for (int j = 0; j < 3; j++) {
      if (std::isnan(V(i, j)))
        has_nan = true;
    }

    if (has_nan) {
      cout << "Vertex " << i << " has nans. " << V.row(i) << endl;
      result = false;
      vertex_issues.insert(i);
    }
  }

  cout << "Minimum face area: " << min_area << endl;
  cout << "Number unused vertices: " << unused << endl;

  return result;
}

bool sparseMatrixHasNaN(const Eigen::SparseMatrix<double> &A) {
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      double x = it.value();

      if (std::isnan(x))
        return true;
    }
  }

  return false;
}

void writeMeshWithMarkers(const char* fn,
                          const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &M) {
  char* fn_off = new char[strlen(fn) + 15];
  sprintf(fn_off, "%s.off", fn);
  igl::writeOFF(fn_off, V, F);
  sprintf(fn_off, "%s_orig.txt", fn);
  FILE* of = fopen(fn_off, "w");
  for (int i = 0; i < M.rows(); ++i)
    fprintf(of, "%d\n", M(i));
  fclose(of);

  delete[] fn_off;
}

} // namespace Helpers
