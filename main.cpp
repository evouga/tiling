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

using namespace std;

// Input/output is the same!!
void removeDuplicates(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &orig) {
  Eigen::MatrixXd Vo;
  Eigen::MatrixXi Fo;
  Eigen::VectorXi I, origo;
  igl::remove_duplicates(V, F, Vo, Fo, I, 1e-6);
  origo.resize(Vo.rows());
  for (int i = 0; i < orig.rows(); ++i) {
    origo(I(i)) = orig(i);
  }
  V = Vo;
  F = Fo;
  orig = origo;
}

// Input/output is the same!!
void improveMesh(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &orig, float triangle_area) {
  Eigen::MatrixXi Fo;
  igl::collapse_small_triangles(V, F, triangle_area, Fo);
  F = Fo;
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
  /*
  for (int i = 0; i < triOrig.rows(); ++i) {
    printf("%d -> %d\n", i, J(i));
    triOrig(i) = orig(J(i));
  }
  */

  return;

  // Comment the rest of this out.
  igl::copyleft::cgal::outer_hull(V,F, Vo,Fo, J,flip);

  origo.resize(Vo.rows());
  origo = Eigen::VectorXi::Constant(Vo.rows(), GLOBAL::nonoriginal_marker);

  // Need to set the original, based off of the original from the former.
  int no = 0;
  for (int i = 0; i < Fo.rows(); ++i) {
    // Get the row from F
    const auto& Frow = F.row(J(i));

    for (int j = 0; j < Fo.cols(); ++j) {
      if (orig(Frow(j)) == GLOBAL::original_marker) {
        origo(Fo(i, j)) = GLOBAL::original_marker;
        no++;
      }
    }
  }

  triV = Vo;
  triF = Fo;
  triOrig = origo;


  fprintf(stderr,"Verts now %d (vs %d), F:%d (vs %d), o:%d\n",
         triV.rows(),Vo.rows(), triF.rows(), Fo.rows(), triOrig.rows());

  /*
  // Remove the bottom
  double zmin = Vo(0,2);
  for (int i = 0; i < Vo.rows(); ++i) {
    if (Vo(i, 2) < zmin) zmin = Vo(i, 2);
  }
  Eigen::VectorXi newIdx(Vo.rows());
  Eigen::VectorXi usedV(Vo.rows());
  usedV.setZero();
  for (int i = 0; i < Fo.rows(); ++i) {
    for (int j = 0; j < Fo.cols(); ++j) {
      // If we're at the min
      if (Vo(Fo(i, j)) == zmin) {
        // Check to see if any of our neighbors are not the min.
        bool valid = false;
        for (int k = 0; k < Fo.cols(); ++k) {
          if (Vo(Fo(i, (j + k) % Fo.cols())) != zmin) {
            valid = true;
            break;
          }
        }

        if (valid) {
          usedV(Fo(i, j)) = 1;
        }
      } else {
        usedV(Fo(i, j)) = 1;
      }
    }
  }

  // Create the new V vector and the map from indices of one to indices of the other.
  int tempV = 0;
  for (int i = 0; i < Vo.rows(); ++i) {
    if (usedV(i)) tempV++;
  }
  triV.resize(tempV, 3);
  // Also update triOrig
  triOrig.resize(tempV);
  tempV = 0;
  Eigen::VectorXi Vmap(Vo.rows());
  for (int i = 0; i < Vo.rows(); ++i) {
    if (usedV(i)) {
      triV.row(tempV) = Vo.row(i);
      triOrig(tempV) = origo(i);
      Vmap(i) = tempV++;
    } else {
      Vmap(i) = -1;
    }
  }
  // Also remove unused faces
  int tempF = 0;
  for (int i = 0; i < Fo.rows(); ++i) {
    bool used = true;
    for (int j = 0; j < Fo.cols(); ++j) {
      if (!usedV(Fo(i, j))) {
        used = false;
        break;
      }
    }
    if (used) {
      tempF++;
    }
  }
  triF.resize(tempF, Fo.cols());
  tempF = 0;
  for (int i = 0; i < Fo.rows(); ++i) {
    bool used = true;
    for (int j = 0; j < Fo.cols(); ++j) {
      if (!usedV(Fo(i, j))) {
        used = false;
        break;
      }
    }
    if (used) {
      triF.row(tempF) = Fo.row(i);
      tempF++;
    }
  }
  // Update the inidices of F.
  for (int i = 0; i < triF.rows(); ++i) {
    for (int j = 0; j < triF.cols(); ++j) {
      // Get the new V
      triF(i, j) = Vmap(triF(i, j));
    }
  }

  printf("Sizes: V:%d(vs %d) F:%d no:%d\n", triV.rows(), Vo.rows(), triF.rows(), no);
  */
}


int nComponents(const Eigen::VectorXi &comps) {
  std::set<int> u;
  for (int i = 0; i < comps.rows(); ++i) {
    u.insert(comps[i]);
  }
  return u.size();
}

void combineMeshes(const Eigen::MatrixXd &bV, const Eigen::MatrixXi &bF,
                   const Eigen::MatrixXd &tV, const Eigen::MatrixXi &tF,
                   Eigen::MatrixXd &V, Eigen::MatrixXi &F, bool shift = true) {

  auto max_bot = bV.colwise().maxCoeff();
  auto min_top = tV.colwise().minCoeff();
  double shift_z = max_bot(2) - min_top(2);
  // Shift the top vertices.
  Eigen::MatrixXd temptV = tV;
  if (shift) {
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
                   const Eigen::MatrixXd &tV, const Eigen::MatrixXi &tF,
                   Eigen::MatrixXd &V, Eigen::MatrixXi &F, bool shift = true) {
  // Combine the two of them.
  V.resize(bV.rows() + tV.rows(), 3);
  F.resize(bF.rows() + tF.rows(), 3);
  auto max_bot = bV.colwise().maxCoeff();
  auto min_top = tV.colwise().minCoeff();
  double shift_z = max_bot(2) - min_top(2);
  printf("Shift is %lf (%lf-%lf)\n", shift_z, max_bot(2), min_top(2));

  V << bV, tV;
  F << bF, tF;
  // Change the top vertices so they have the appropriate z-values.
  for (int i = bV.rows(); i < V.rows(); ++i) {
    if (shift) V(i, 2) += shift_z;
  }
  // Change the face indices so they have the correct vertex idx.
  for (int i = bF.rows(); i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      // Increment it to point to the correct location.
      F(i, j) += bV.rows();
    }
  }

	Eigen::MatrixXd Vu;
  Eigen::MatrixXi Fu;
  Eigen::VectorXi I;
  igl::remove_duplicates(V,F, Vu,Fu, I);

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
void getOffsetSurface(int bot_slice_no, SliceStack &ss,
                      Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces,
                      Eigen::VectorXi &botorig,
                      Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
                      Eigen::VectorXi &toporig,
                      Eigen::MatrixXd &offsetV, Eigen::MatrixXi &offsetF,
                      Eigen::VectorXi &orig,
                      bool view=false) {
  igl::viewer::Viewer v;

  ss.triangulateSlice(bot_slice_no, 0.005,
                      botverts, botfaces, topverts, topfaces,
											botorig, toporig);

  if (view) {
    Eigen::MatrixXd C;
    igl::jet(botorig, true, C);
    v.data.set_mesh(botverts, botfaces);
    v.data.set_face_based(true);
    v.data.set_colors(C);
    for (int i = 0; i < botorig.rows(); ++i) {
      if (botorig(i) == 0) {
        v.data.add_label(botverts.row(i), "o");
      } else if (botorig(i) == 1) {
        v.data.add_label(botverts.row(i), "x");
      } else {
        v.data.add_label(botverts.row(i), to_string(botorig(i)));
      }
    }
    Eigen::MatrixXd allverts;
    Eigen::MatrixXi allfaces;
    combineMeshes_2(botverts, botfaces, topverts, topfaces,
                  allverts, allfaces, false);
    v.data.clear();
    v.data.set_mesh(allverts, allfaces);
    v.data.set_face_based(true);
    v.launch();
  }

  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT, TF;
	Eigen::VectorXi TO;
  ss.tetrahedralizeSlice(botverts, botfaces, topverts, topfaces,
												 botorig, toporig,
												 TV, TT, TF, TO);

  cout << "Triangulated tile contains ";
  cout << botverts.rows() << " verts on bottom face and ";
  cout << topverts.rows() << " verts on top face" << endl;

  Eigen::VectorXd Z;
	ss.computeLaplace(bot_slice_no, TV, TT, TF, TO, Z);
  if (view) TetMeshViewer::viewTetMesh(TV, TT, TF, Z, true);

  //if (view) TetMeshViewer::viewOffsetSurface(TV, TF, TT, Z);

  // Now, get the offset surface
  double offset = 0.5;
  //OffsetSurface::Triangulation T;
  //generateOffsetSurface(TV, TT, Z, offset, offsetV, offsetF, T);
  OffsetSurface::generateOffsetSurface_naive(
      TV, TT, Z, offset, offsetV, offsetF);
  if (view) {
    v.data.set_mesh(offsetV, offsetF);
    v.launch();
  }

  Eigen::VectorXi comp;
  igl::components(offsetF, comp);
  int nc = nComponents(comp);
  while (nc > 1) {
    printf("With an offset of %lf, number of components is %d\n",
           offset, nc);
    offset = (offset + 1.0) / 2.0;
    //generateOffsetSurface(T, offset, offsetV, offsetF);
    OffsetSurface::generateOffsetSurface_naive(
        TV, TT, Z, offset, offsetV, offsetF);
    igl::components(offsetF, comp);
    nc = nComponents(comp);
    if (view) {
      v.data.clear();
      v.data.set_mesh(offsetV, offsetF);
      v.data.set_face_based(true);
      v.launch();
    }
  }

  auto maxs = offsetV.colwise().maxCoeff();
  auto mins = offsetV.colwise().minCoeff();
  orig.resize(offsetV.rows());
  // It's original if it's on the boundary.
  for (int i = 0; i < offsetV.rows(); ++i) {
    if (offsetV(i,2) == maxs(2) ||
        offsetV(i,2) == mins(2)) {
      orig(i) = GLOBAL::original_marker;
    } else {
      orig(i) = GLOBAL::nonoriginal_marker;
    }
  }
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    cerr << "Must specify the base filename and object name" << endl;
    return -1;
  }

  SliceStack ss(argv[1], argv[2]);

  cout << "Loaded " << ss.getNumSlices() << " slices" << endl;
  int good_start = -1;
  for (int i = 0; i < ss.getNumSlices(); ++i) {
    if (ss.getSizeAt(i) > 0) {
      good_start = i;
      printf("Using slice number %d\n", good_start);
      break;
    }
  }
  if (good_start < 0) {
    printf("ERROR: Couldn't find valid slices!\n");
    return -1;
  }
  Eigen::MatrixXd bV, tV;
  Eigen::MatrixXi bF, tF;
  Eigen::VectorXi borig, torig;

  Eigen::MatrixXd topverts, botverts;
  Eigen::MatrixXi topfaces, botfaces;
  Eigen::VectorXi toporig, botorig;

  getOffsetSurface(good_start, ss,
                   botverts, botfaces, botorig,
                   topverts, topfaces, toporig,
                   bV, bF, borig, true);
  good_start++;
  // Next time's bottom will be last time's top.
  botverts = topverts;
  botfaces = topfaces;
  botorig = toporig;

  // Get the next cube, using stuff from last time.
  getOffsetSurface(good_start, ss,
                   botverts, botfaces, botorig,
                   topverts, topfaces, toporig,
                   tV, tF, torig, true);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  combineMeshes(bV, bF, tV, tF, V, F);
  Eigen::VectorXi orig(V.rows());
  auto mins = bV.colwise().minCoeff();
  auto mids = bV.colwise().maxCoeff();
  auto maxs = V.colwise().maxCoeff();

  int n_orig = 0;
  for (int i = 0; i < V.rows(); ++i) {
    if (V(i, 2) == mins(2) ||
        V(i, 2) == mids(2) ||
        V(i, 2) == maxs(2)) {
      orig(i) = GLOBAL::original_marker;
      ++n_orig;
    } else {
      orig(i) = GLOBAL::nonoriginal_marker;
    }
  }
  printf("Number orig: %d/%ld\n", n_orig, V.rows());

  // Display them
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.launch();

  // Let's write it off right now.
  igl::writeOFF("both_slices_pre.off", V, F);
  FILE* ofpre = fopen("both_slices_pre_orig.txt", "w");
  for (int i = 0; i < orig.rows(); ++i) {
    fprintf(ofpre, "%d\n", orig(i));
  }
  fclose(ofpre);

  Eigen::MatrixXd Vborder;
  Eigen::MatrixXi Fborder;
  Eigen::VectorXi Oborder;
  extractShell(V, F, orig, Vborder, Fborder, Oborder);
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();
  fprintf(stderr,"extracted shell\n");
  removeDuplicates(Vborder, Fborder, Oborder);
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();
  fprintf(stderr, "removed duplicates\n");
  improveMesh(Vborder, Fborder, Oborder, 0.01);
  viewer.data.clear();
  viewer.data.set_mesh(Vborder, Fborder);
  viewer.launch();
  fprintf(stderr, "improved mesh\n");
  igl::writeOFF("both_slices.off", Vborder, Fborder);
  // Also write the original vertices
  FILE* of = fopen("both_slices_orig.txt", "w");
  for (int i = 0; i < Oborder.rows(); ++i) {
    fprintf(of, "%d\n", Oborder(i));
  }
  fclose(of);

  printf("Written to output.\n");
  Eigen::MatrixXd Vcurve;
  //computeCurvatureFlow(V, F, orig, 0.1, Vcurve);
  //viewer.data.clear();
  //viewer.data.set_mesh(Vcurve, F);
  //viewer.launch();
  //biharmonic(V, F, orig, 0.1, Vcurve);
  biharmonic(Vborder, Fborder, Oborder, 0.1, Vcurve);

	/*
  Eigen::MatrixXd V(botverts.rows() + topverts.rows(), 3);
  for(int i=0; i<botverts.rows(); i++)
    V.row(i) = botverts.row(i);
  for(int i=0; i<topverts.rows(); i++)
    V.row(botverts.rows() + i) = topverts.row(i);
  Eigen::MatrixXi F(botfaces.rows() + topfaces.rows(), 3);
  for(int i=0; i<botfaces.rows(); i++)
    F.row(i) = botfaces.row(i);
  for(int i=0; i<topfaces.rows(); i++)
  {
    for(int j=0; j<3; j++)
      F(botfaces.rows() + i, j) = topfaces(i,j) + botverts.rows();
  }

  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_face_based(true);
  viewer.launch();
	*/
}
