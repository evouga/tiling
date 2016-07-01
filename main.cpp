#include <set>

#include <igl/components.h>
#include <igl/jet.h>
#include <igl/unique.h>
#include <igl/viewer/Viewer.h>

#include "curvatureFlow.h"
#include "glob_defs.h"
#include "offsetSurface.h"
#include "SliceStack.h"
#include "viewTetMesh.h"

using namespace std;

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
  // Combine the two of them.
  V.resize(bV.rows() + tV.rows(), 3);
  F.resize(bF.rows() + tF.rows(), 3);
  auto max_bot = bV.colwise().maxCoeff();
  auto min_top = tV.colwise().minCoeff();
  double shift_z = max_bot(2) - min_top(2);
  printf("Shift is %lf (%lf-%lf)\n", shift_z, max_bot(2), min_top(2));

  V << bV, tV;
  F << bF, tF;
  // Change the vertices so they have the appropriate z-values.
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
}

void getOffsetSurface(int bot_slice_no, SliceStack &ss,
                      Eigen::MatrixXd &offsetV, Eigen::MatrixXi &offsetF,
                      Eigen::VectorXi &orig,
                      bool view=false) {

  igl::viewer::Viewer v;

  Eigen::MatrixXd botverts;
  Eigen::MatrixXd topverts;
  Eigen::MatrixXi botfaces;
  Eigen::MatrixXi topfaces;
	Eigen::VectorXi toporig;
	Eigen::VectorXi botorig;
  //                                      don't scale
  ss.triangulateSlice(bot_slice_no, 0.05, false,
                      botverts, botfaces, topverts, topfaces,
											botorig, toporig);

  if (view) {
    Eigen::MatrixXd C;
    igl::jet(toporig, true, C);
    v.data.set_mesh(topverts, topfaces);
    v.data.set_face_based(true);
    //v.data.set_colors(C);
    for (int i = 0; i < toporig.rows(); ++i) {
      if (toporig(i) == 0) {
        v.data.add_label(topverts.row(i), "o");
      } else if (toporig(i) == 1) {
        v.data.add_label(topverts.row(i), "x");
      } else {
        v.data.add_label(topverts.row(i), to_string(toporig(i)));
      }
    }
    v.launch();

    Eigen::MatrixXd allverts;
    Eigen::MatrixXi allfaces;
    combineMeshes(botverts, botfaces, topverts, topfaces,
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

  cout << "Triangulated tile contains " << botverts.rows() << " verts on bottom face and " << topverts.rows() << " verts on top face" << endl;  

  Eigen::VectorXd Z;
	ss.computeLaplace(bot_slice_no, TV, TT, TF, TO, Z);
  if (view) TetMeshViewer::viewTetMesh(TV, TT, TF, Z, true);

  //if (view) TetMeshViewer::viewOffsetSurface(TV, TF, TT, Z);

  // Now, get the offset surface
  OffsetSurface::Triangulation T;
  double offset = 0.5;
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

int main(int argc, char *argv[])
{
  if(argc != 3)
  {
    cerr << "Must specify the base filename and object name" << endl;
    return -1;
  }

  SliceStack ss(argv[1], argv[2]);

  cout << "Loaded " << ss.getNumSlices() << " slices" << endl;
  int good_start = 0;
  for (int i = 0; i < ss.getNumSlices(); ++i) {
    if (ss.getSizeAt(i) > 0) {
      good_start = i;
      break;
    }
  }
  printf("Using slice number %d\n", good_start);
  Eigen::MatrixXd bV, tV;
  Eigen::MatrixXi bF, tF;
  Eigen::VectorXi borig, torig;
  getOffsetSurface(good_start, ss, bV, bF, borig, false);
  good_start++;
  getOffsetSurface(good_start, ss, tV, tF, torig, false);

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
  printf("Number orig: %d/%d\n", n_orig, V.rows());

  // Display them
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.launch();

  Eigen::MatrixXd Vcurve;
  computeCurvatureFlow(V, F, orig, 0.1, Vcurve);
  viewer.data.clear();
  viewer.data.set_mesh(Vcurve, F);
  viewer.launch();
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
