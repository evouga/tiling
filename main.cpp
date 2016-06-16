#include <igl/viewer/Viewer.h>
#include <igl/jet.h>

#include "SliceStack.h"
#include "viewTetMesh.h"

using namespace std;

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

  Eigen::MatrixXd botverts;
  Eigen::MatrixXd topverts;
  Eigen::MatrixXi botfaces;
  Eigen::MatrixXi topfaces;
	Eigen::VectorXi toporig;
	Eigen::VectorXi botorig;
  ss.triangulateSlice(good_start, 0.05, botverts, botfaces, topverts, topfaces,
											botorig, toporig);

  /*
  igl::viewer::Viewer viewer;
	Eigen::MatrixXd C;
	igl::jet(toporig, true, C);
  viewer.data.set_mesh(topverts, topfaces);
  viewer.data.set_face_based(true);
	//viewer.data.set_colors(C);
	for (int i = 0; i < toporig.rows(); ++i) {
		if (toporig(i) == 0) {
			viewer.data.add_label(topverts.row(i), "o");
		} else if (toporig(i) == 1) {
			viewer.data.add_label(topverts.row(i), "x");
		} else {
			viewer.data.add_label(topverts.row(i), to_string(toporig(i)));
		}
	}
  viewer.launch();

	igl::jet(botorig, true, C);
	viewer.data.clear();
  viewer.data.set_mesh(botverts, botfaces);
  viewer.data.set_face_based(true);
	viewer.data.set_colors(C);
  viewer.launch();
  */

  Eigen::MatrixXd TV;
  Eigen::MatrixXi TT, TF;
	Eigen::VectorXi TO;
  ss.tetrahedralizeSlice(botverts, botfaces, topverts, topfaces,
												 botorig, toporig,
												 TV, TT, TF, TO);

  TetMeshViewer::viewTetMesh(TV, TT, TF);
  cout << "Triangulated tile contains " << botverts.rows() << " verts on bottom face and " << topverts.rows() << " verts on top face" << endl;  

	//ss.computeLaplace(TV, TT, TF, TO);
	ss.computeLaplace(good_start, TV, TT, TF, TO);
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
