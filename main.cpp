#include <set>

#include <igl/boundary_facets.h>
#include <igl/components.h>
#include <igl/collapse_small_triangles.h>
#include <igl/copyleft/cgal/outer_hull.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/jet.h>
#include <igl/remove_duplicates.h>
#include <igl/remove_unreferenced.h>
#include <igl/resolve_duplicated_faces.h>
#include <igl/simplify_polyhedron.h>
#include <igl/unique.h>
#include <igl/viewer/Viewer.h>
#include <igl/writeOFF.h>

#include "curvatureFlow.h"
#include "glob_defs.h"
#include "offsetSurface.h"
#include "SliceStack.h"
#include "viewTetMesh.h"
#include "Helpers.h"

using namespace std;

// Input/output is the same!!
void improveMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                 Eigen::VectorXi &orig, float triangle_area=0.01) {
  Eigen::MatrixXi Fo;
  igl::collapse_small_triangles(V, F, triangle_area, Fo);
  F = Fo;
}

// Returns the genus of the given mesh. Genus can be calculated with Euler's formula:
//   2 (1 - g) = #V - #E + #F
// 
// Given are #V and #F, all that is needed is computing #E
int getGenus(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
  // The number of edges is (F.rows() * 3) / 2
  return (V.rows() - 0.5 * F.rows() - 2) / 2;
}

void extractShell(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &orig,
                  Eigen::MatrixXd &triV, Eigen::MatrixXi &triF, Eigen::VectorXi &triOrig) {
  Eigen::MatrixXd Vo;
  Eigen::MatrixXi Fo;
  Eigen::VectorXi origo;

  Eigen::VectorXi J, flip;
  igl::resolve_duplicated_faces(F, Fo, J);
  // J is the same size as V
  igl::remove_unreferenced(V, Fo, triV, triF, J);
  triOrig.resize(triV.rows());
  triOrig.setZero();
  for (int i = 0; i < orig.rows(); ++i) {
    //printf("%d -> %d\n", i, J(i));
    if (J(i) != -1) {
      triOrig(J(i)) = orig(i);
    }
  }
}


int nComponents(const Eigen::VectorXi &comps) {
  std::set<int> u;
  for (int i = 0; i < comps.rows(); ++i) {
    u.insert(comps[i]);
  }
  return u.size();
}

void combineMeshes(const Eigen::MatrixXd &bV, const Eigen::MatrixXi &bF,
                   const Eigen::VectorXi &bO,
                   const Eigen::MatrixXd &tV, const Eigen::MatrixXi &tF,
                   const Eigen::VectorXi &tO,
    Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O,
    bool shift = true) {

  Eigen::MatrixXd temptV = tV;
  if (shift) {
    // Put the two of them together (so that the bax of the bottom and the min
    // of the top are the same
    auto max_bot = bV.colwise().maxCoeff();
    auto min_top = tV.colwise().minCoeff();
    double shift_z = max_bot(2) - min_top(2);
    // Shift the top vertices.
    for (int i = 0; i < temptV.rows(); ++i) {
      temptV(i, 2) += shift_z;
    }
  }

  // Then use CGAL's MESH_BOOLEAN function.
  Eigen::VectorXi J;
  igl::copyleft::cgal::mesh_boolean(bV,bF, temptV,tF,
                                    igl::MeshBooleanType::MESH_BOOLEAN_TYPE_UNION,
                                    V,F, J);
}

void combineMeshes_2(const Eigen::MatrixXd &bV, const Eigen::MatrixXi &bF,
                     const Eigen::VectorXi &bO,
                     const Eigen::MatrixXd &tV, const Eigen::MatrixXi &tF,
                     const Eigen::VectorXi &tO,
                     Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O) {
  // Combine the two of them.
  V.resize(bV.rows() + tV.rows(), 3);
  F.resize(bF.rows() + tF.rows(), 3);
  O.resize(bO.rows() + tO.rows());

  V << bV, tV;
  F << bF, tF;
  O << bO, tO;

  // Change the face indices so they have the correct vertex idx.
  for (int i = bF.rows(); i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      F(i, j) += bV.rows();
    }
  }

  Helpers::removeDuplicates(V, F, O);

  /*
  // Make sure the vertices are all unique
	// Get all the unique vertices
	Eigen::MatrixXd Vu;
	// contains mapping from unique to all indices
	// Size: #V
	Eigen::VectorXi unique_to_all;
	// contains mapping from all indices to unique
	// Size: #V_rep
	Eigen::VectorXi all_to_unique;
	igl::unique_rows(V, Vu,unique_to_all,all_to_unique);
  // Remember these.
  V = Vu;

  // Also need to update faces.
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      F(i, j) = all_to_unique(F(i, j));
    }
  }
  // And remove duplicate faces.
  Eigen::MatrixXi Fu;
  igl::unique_rows(F, Fu,unique_to_all,all_to_unique);
  F = Fu;
  */
}

// Can pass in topverts, topfaces, and toporig if you want to use them.
void getOffsetSurface(
    int bot_slice_no, SliceStack &ss,
    Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, Eigen::VectorXi &botorig,
    Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces, Eigen::VectorXi &toporig,
    Eigen::MatrixXd &offsetV, Eigen::MatrixXi &offsetF, Eigen::VectorXi &orig,
    const std::vector<int> &allowed_bot, const std::vector<int> &allowed_top,
    bool view=false) {
  igl::viewer::Viewer v;

  ss.triangulateSlice(bot_slice_no, GLOBAL::TRI_AREA,
                      botverts, botfaces, topverts, topfaces,
											botorig, toporig, allowed_bot, allowed_top);

  if (view) {
    Eigen::MatrixXd C;
    igl::jet(botorig, true, C);
    v.data.clear();
    v.data.set_mesh(botverts, botfaces);
    v.data.set_face_based(true);
    v.data.set_colors(C);
    for (int i = 0; i < botorig.rows(); ++i) {
      if (botorig(i) == GLOBAL::nonoriginal_marker) {
        v.data.add_label(botverts.row(i), "o");
      } else if (botorig(i) == GLOBAL::original_marker) {
        v.data.add_label(botverts.row(i), "x");
      } else {
        v.data.add_label(botverts.row(i), to_string(botorig(i)));
      }
    }
    Eigen::MatrixXd allverts;
    Eigen::MatrixXi allfaces;
    Eigen::VectorXi allorig;
    combineMeshes_2(botverts, botfaces, botorig,
                    topverts, topfaces, toporig,
                    allverts, allfaces, allorig);
    v.data.clear();
    v.data.set_mesh(allverts, allfaces);
    v.data.set_face_based(true);
    v.launch();
  }

  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT, TF;
	Eigen::VectorXi TO;
  ss.tetrahedralizeSlice(botverts, botfaces, topverts, topfaces,
												 botorig, toporig, TV, TT, TF, TO);

  cout << "Triangulated tile contains ";
  cout << botverts.rows() << " verts on bottom face and ";
  cout << topverts.rows() << " verts on top face" << endl;

  Eigen::VectorXd Z;
	ss.computeLaplace(TV, TT, TF, TO, Z);
  if (view) TetMeshViewer::viewTetMesh(TV, TT, TF, Z, true);

  //if (view) TetMeshViewer::viewOffsetSurface(TV, TF, TT, Z);

  // Now, get the offset surface
  double offset = 0.5;
  //OffsetSurface::Triangulation T;
  //generateOffsetSurface(TV, TT, Z, offset, offsetV, offsetF, T);
  OffsetSurface::generateOffsetSurface_naive(
      TV, TT, TO, Z, offset, offsetV, offsetF, orig);
  if (view) {
    v.data.clear();
    v.data.set_mesh(offsetV, offsetF);
    Eigen::MatrixXd origCols;
    igl::jet(orig, true, origCols);
    v.data.set_colors(origCols);
    v.launch();
  }

  Eigen::VectorXi comp;
  igl::components(offsetF, comp);
  int nc = nComponents(comp);
  while (nc > 1) {
    printf("With an offset of %lf, number of components is %d\n", offset, nc);
    offset = (offset + 1.0) / 2.0;
    //generateOffsetSurface(T, offset, offsetV, offsetF);
    OffsetSurface::generateOffsetSurface_naive(
        TV, TT, TO, Z, offset, offsetV, offsetF, orig);
    igl::components(offsetF, comp);
    nc = nComponents(comp);
    if (view) {
      v.data.clear();
      v.data.set_mesh(offsetV, offsetF);
      v.data.set_face_based(true);
      Eigen::MatrixXd origCols;
      igl::jet(orig, true, origCols);
      v.data.set_colors(origCols);
      v.launch();
    }
  }

  if (view) {
    // Some debug stuff
    igl::writeOFF("debug_tets.off", TV, TF);
    // Now need to write the tets and heats as well
    FILE* tet_of = fopen("debug_tets.txt", "w");
    for (int i = 0; i < TT.rows(); ++i) {
      for (int j = 0; j < TT.cols(); ++j) {
        fprintf(tet_of, "%d ", TT(i, j));
      }
      fprintf(tet_of, "\n");
    }
    fclose(tet_of);
    // Heats
    FILE* heat_of = fopen("debug_heats.txt", "w");
    for (int i = 0; i < Z.rows(); ++i) {
      fprintf(heat_of, "%lf\n", Z(i));
    }
    fclose(heat_of);
    // Origs
    FILE* orig_of = fopen("debug_orig.txt", "w");
    for (int i = 0; i < TO.rows(); ++i) {
      fprintf(orig_of, "%d\n", TO(i));
    }
    fclose(orig_of);
  }

  /*
  auto maxs = offsetV.colwise().maxCoeff();
  auto mins = offsetV.colwise().minCoeff();
  orig.resize(offsetV.rows());
  // It's original if it's on the boundary.
  for (int i = 0; i < offsetV.rows(); ++i) {
    if (offsetV(i,2) == maxs(2) || offsetV(i,2) == mins(2)) {
      orig(i) = GLOBAL::original_marker;
    } else {
      orig(i) = GLOBAL::nonoriginal_marker;
    }
  }
  */
}

int main(int argc, char *argv[]) {
  if(argc < 3) {
    fprintf(stderr, "usage: %s <BASE_FILENAME> <CONTOUR_NAME> [<good_start> <num_slices>]\n",
            argv[0]);
    return -1;
  }

  SliceStack ss(argv[1], argv[2]);

  cout << "Loaded " << ss.getNumSlices() << " slices" << endl;
  int good_start = 0;
  int num_slices = 10;
  if (argc == 5) {
    good_start = atoi(argv[3]);
    num_slices = atoi(argv[4]);
  }

  while (ss.getSizeAt(good_start) == 0) {
    good_start++;
    if (good_start > ss.getNumSlices()) {
      // Quit now.
      printf("ERROR: Couldn't find valid slices!\n");
      return -1;
    }
  }
  printf("Starting with slice number %d\n", good_start);

  // Use an empty vector to allow everything.
  //std::vector<int> all_allowed = {0};
  std::vector<int> all_allowed;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi orig;
  Eigen::MatrixXd bV, tV;
  Eigen::MatrixXi bF, tF;
  Eigen::VectorXi borig, torig;

  Eigen::MatrixXd topverts, botverts;
  Eigen::MatrixXi topfaces, botfaces;
  Eigen::VectorXi toporig, botorig;
  getOffsetSurface(good_start, ss,
                   botverts, botfaces, botorig,
                   topverts, topfaces, toporig,
                   bV, bF, borig, all_allowed, all_allowed, true);
  // Make sure this is consistent if they only ask for a single slice.
  V = bV;
  F = bF;
  orig = borig;

  good_start++;
  int size = ss.getSizeAt(good_start);
  int num = 0;
  while (size > 0 && ++num < num_slices) {
    // Next time's bottom will be last time's top.
    botverts = topverts;
    botfaces = topfaces;
    botorig = toporig;

    // Get the next cube, using stuff from last time.
    getOffsetSurface(good_start, ss,
                     botverts, botfaces, botorig,
                     topverts, topfaces, toporig,
                     tV, tF, torig, all_allowed, all_allowed, false /* display */);


    combineMeshes_2(bV, bF, borig, tV, tF, torig, V, F, orig);
    bV = V;
    bF = F;
    borig = orig;
  
    // Get the size of the next contour.
    good_start++;
    size = ss.getSizeAt(good_start);
  }

  // Display them
  igl::viewer::Viewer viewer;
  Eigen::MatrixXd origCols;
  igl::jet(orig, true, origCols);
  viewer.data.clear();
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(origCols);
  viewer.launch();

  Eigen::MatrixXd Vborder;
  Eigen::MatrixXi Fborder;
  Eigen::VectorXi Oborder;
  if (num > 1) {
    fprintf(stderr, "Extracting shell...\n");
    extractShell(V, F, orig, Vborder, Fborder, Oborder);
  } else {
    // If we only use a single thing, we need to take a little special extra care.
    Vborder = V;
    Fborder = F;
    Oborder = orig;
    extractShell(V, F, orig, Vborder, Fborder, Oborder);

    fprintf(stderr, "Removing duplicate vertices (%lu verts)\n", Vborder.rows());
    Helpers::removeDuplicates(Vborder, Fborder, Oborder);
  }
  /*
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();
  fprintf(stderr,"extracted shell\nRemoving duplicates...\n");
  removeDuplicates(Vborder, Fborder, Oborder);
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();

  // This isn't working... But it's also not (really) necessary.
  fprintf(stderr, "removed duplicates\nimproving mesh...");
  improveMesh(Vborder, Fborder, Oborder);
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();
  fprintf(stderr, "improved mesh\n");
  */
  igl::writeOFF("both_slices.off", Vborder, Fborder);
  // Also write the original vertices
  FILE* of = fopen("both_slices_orig.txt", "w");
  for (int i = 0; i < Oborder.rows(); ++i) {
    fprintf(of, "%d\n", Oborder(i));
  }
  fclose(of);

  printf("Written to output [both_slices.off] and [both_slices_orig.txt].\n");
  Eigen::MatrixXd Vcurve;
  //computeCurvatureFlow(V, F, orig, 0.1, Vcurve);
  //viewer.data.clear();
  //viewer.data.set_mesh(Vcurve, F);
  //viewer.launch();
  //biharmonic(V, F, orig, 0.1, Vcurve);
  biharmonic_view(Vborder, Fborder, Oborder, Vcurve, false);
}
