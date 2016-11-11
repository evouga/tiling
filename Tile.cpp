#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>

#include <igl/jet.h>
#include <igl/triangle/triangulate.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "Tile.h"
#include "Slice.h"

using namespace Eigen;
using namespace std;

#define ZSCALE_HACK

const double Tile::tilePadding_ = 0.10;

Tile::Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox) 
    : bottom_(bottom), top_(top), bbox_(bbox)
{ }

Tile::~Tile()
{
}

/*
// Will resize V and E to fit
void Tile::getOrig(Eigen::MatrixXd &Vtop, Eigen::MatrixXd &Vbot) {
	unsigned int total_pts = bottom_.getNumPts() + top_.getNumPts();
	Vtop.resize(top_.getNumPts(), 3);
	Vbot.resize(bottom_.getNumPts(), 3);
	Eigen::MatrixXi E;
	E.resize(bottom_.getNumPts(), 3);
	int offset = addOrig(bottom_, Vbot, E, 0);
	assert(offset == bottom_.getNumPts());
	fprintf(stderr, "offset is %d and size is %d v is %d\n",
					offset, bottom_.getNumPts(), Vbot.rows());
	// Add the z-index
	for (int i = 0; i < offset; ++i) {
		Vbot(i, 2) = -GLOBAL::z_lim;
	}

	E.resize(top_.getNumPts(), 3);
	offset = addOrig(top_, Vtop, E, 0);
	for (int i = 0; i < offset; ++i) {
		Vtop(i, 2) = GLOBAL::z_lim;
	}
}
*/

// Adds original points to V and E. Should be already the correct size
int Tile::addOrig(const Slice &s, 
									Eigen::MatrixXd &V, Eigen::MatrixXi &E,
									Eigen::VectorXi &VM, Eigen::VectorXi &EM,
									Eigen::MatrixXd &lims, int offset) {
  bool first = true;
  lims.resize(2, 2);
	int numcontours = s.contours.size();
	for(int i=0; i<numcontours; i++)
	{
		int numpts = s.contours[i].x.size();
		for(int j=0; j<numpts; j++)
		{
			Vector2d pt(s.contours[i].x[j], s.contours[i].y[j]);
			V(offset+j, 0) = pt(0);
			V(offset+j, 1) = pt(1);
			int prev = (j == 0 ? numpts-1 : j-1);
			E(offset+j, 0) = offset + j;
			E(offset+j, 1) = offset + prev;
			VM(offset+j) = GLOBAL::original_marker;
			EM(offset+j) = GLOBAL::original_marker;
      if (first) {
        lims(0, 0) = lims(0, 1) = pt(0);
        lims(1, 0) = lims(1, 1) = pt(1);
        first = false;
      } else {
        lims(0, 0) = min(lims(0, 0), pt(0));
        lims(0, 1) = max(lims(0, 1), pt(0));
        lims(1, 0) = min(lims(1, 0), pt(1));
        lims(1, 1) = max(lims(1, 1), pt(1));
      }
		}
		offset += numpts;
	}
	return offset;
}

void flood_fill(const Eigen::MatrixXi &faces, Eigen::VectorXi &orig) {
	bool dirty = true;
	while(dirty) {
		dirty = false;
		for (int i = 0; i < faces.rows(); ++i) {
			if (orig(faces(i,0)) == 1 || orig(faces(i,1)) == 1 || orig(faces(i,2)) == 1) {
				for (int j = 0; j < 3; ++j) {
					int vert = faces(i,j);
					if (orig(vert) == 0) {
						orig(vert) = 1;
						dirty = true;
					}
				}
			}
		}
	}

	for (int i = 0; i < orig.rows(); ++i) {
		if (orig(i) == 0) {
			orig(i) = GLOBAL::original_marker;
		}
	}
}

// Triangulates a single contour.
void Tile::triangulateSlice(const Slice &s, 
                            double xmin, double xmax,
                            double ymin, double ymax,
                            double z, double areaBound,
														Eigen::MatrixXd &verts, Eigen::MatrixXi &faces,
														Eigen::VectorXi &orig)
{
	int numcontours = s.contours.size();
	cout << "num contours: " << numcontours << endl;
	int totpts = s.getNumPts();

	MatrixXd V(totpts+4, 2);
	MatrixXi E(totpts+4, 2);
	MatrixXd H(0, 2);
	VectorXi VM(totpts+4);
	VectorXi EM(totpts+4);
  
  /*
  V(0,0) = xmin;
  V(0,1) = ymin;
  V(1,0) = xmax;
  V(1,1) = ymin;
  V(2,0) = xmax;
  V(2,1) = ymax;
  V(3,0) = xmin;
  V(3,1) = ymax;
	E(0,0) = 0;
	E(0,1) = 1;
	E(1,0) = 1;
	E(1,1) = 2;
	E(2,0) = 2;
	E(2,1) = 3;
	E(3,0) = 3;
	E(3,1) = 0;
	VM(0) = GLOBAL::nonoriginal_marker;
	VM(1) = GLOBAL::nonoriginal_marker;
	VM(2) = GLOBAL::nonoriginal_marker;
	VM(3) = GLOBAL::nonoriginal_marker;
	EM(0) = GLOBAL::nonoriginal_marker;
	EM(1) = GLOBAL::nonoriginal_marker;
	EM(2) = GLOBAL::nonoriginal_marker;
	EM(3) = GLOBAL::nonoriginal_marker;
	addOrig(s, V,E,VM,EM, 4);
  */

  Eigen::MatrixXd lims;
	addOrig(s, V,E,VM,EM, lims, 0);
  cout << "Lims are  min:" << endl << V.col(0) << "\nmax:" << endl << V.col(1) << endl;
  cout << "\nvs x " << xmin << "," << xmax;
  cout << " y " << ymin << "," << ymax;
  // Give a gap of 5% of the total area around each vertex.
  auto mins = V.colwise().minCoeff();
  auto maxs = V.colwise().maxCoeff();
  auto gaps = (maxs - mins) * GLOBAL::padding_perc;
  // Add the bounding points after adding the original ones.
  V(totpts,0) = lims(0, 0) - gaps(0);
  V(totpts,1) = lims(1, 0) - gaps(1);
  V(totpts + 1,0) = lims(0, 1) + gaps(0);
  V(totpts + 1,1) = lims(1, 0) - gaps(1);
  V(totpts + 2,0) = lims(0, 1) + gaps(0);
  V(totpts + 2,1) = lims(1, 1) + gaps(1);
  V(totpts + 3,0) = lims(0, 0) - gaps(0);
  V(totpts + 3,1) = lims(1, 1) + gaps(1);
	E(totpts + 0,0) = totpts + 0;
	E(totpts + 0,1) = totpts + 1;
	E(totpts + 1,0) = totpts + 1;
	E(totpts + 1,1) = totpts + 2;
	E(totpts + 2,0) = totpts + 2;
	E(totpts + 2,1) = totpts + 3;
	E(totpts + 3,0) = totpts + 3;
	E(totpts + 3,1) = totpts + 0;
	VM(totpts + 0) = GLOBAL::nonoriginal_marker;
	VM(totpts + 1) = GLOBAL::nonoriginal_marker;
	VM(totpts + 2) = GLOBAL::nonoriginal_marker;
	VM(totpts + 3) = GLOBAL::nonoriginal_marker;
	EM(totpts + 0) = GLOBAL::nonoriginal_marker;
	EM(totpts + 1) = GLOBAL::nonoriginal_marker;
	EM(totpts + 2) = GLOBAL::nonoriginal_marker;
	EM(totpts + 3) = GLOBAL::nonoriginal_marker;
	
  cout << "Rendering triangle here...\n";
	/*
	igl::viewer::Viewer viewer;
	vector<string> labels;
	//Eigen::MatrixXd C;
	//igl::jet(VM, true, C);
	viewer.data.set_mesh(V, E);
  //viewer.data.set_face_based(false);
	//viewer.data.set_colors(C);
	for (int i = 0; i < VM.rows(); ++i) {

		if (VM(i) == 0) {
			viewer.data.add_label(V.row(i), "o");
		} else {
			viewer.data.add_label(V.row(i), "x");
		}
	}
  viewer.launch();
	*/

	stringstream ss;
	ss << "Da" << areaBound << "q";
  
	MatrixXd V2;
	cout << "Delaunay triangulating mesh with " << V.rows() << " verts" << endl;
	igl::triangle::triangulate(V,E,H,VM,EM,ss.str().c_str(),V2,faces,orig);
	cout << "done" << endl;
	verts.resize(V2.rows(), 3);
	flood_fill(faces, orig);
	for(int i=0; i<verts.rows(); i++)
	{
		verts(i, 0) = V2(i,0);
		verts(i, 1) = V2(i,1);
		verts(i, 2) = z;
	}
}

void Tile::triangulateSlices(double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
		Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig)
{
  double xmin = bbox_(0, 0), xmax = bbox_(0, 1),
         ymin = bbox_(1, 0), ymax = bbox_(1, 1);

  double thickness = bottom_.thickness;
  printf("Bounding box is x:%lf,%lf, y:%lf,%lf, z:%lf\n",
         xmin,xmax, ymin,ymax, thickness);
  double xgap = (xmax - xmin) * GLOBAL::padding_perc;
  double ygap = (ymax - ymin) * GLOBAL::padding_perc;
  // If the bottom vertices is not empty, just use what you have.
  if (botverts.rows() == 0) {
    triangulateSlice(bottom_, xmin - xgap, xmax + xgap,
                     ymin - ygap, ymax + ygap, 0,
                     areaBound, botverts, botfaces, bot_orig);
  }
  double prev_z = botverts(0, 2);
#ifdef ZSCALE_HACK
  thickness = (botverts.colwise().maxCoeff() - botverts.colwise().minCoeff()).maxCoeff();
  printf("[%s:%d] Hack in place; thickness is set to %lf, instead of %lf\n",
         __FILE__, __LINE__, thickness, bottom_.thickness);
#endif
  triangulateSlice(top_, xmin - xgap, xmax + xgap,
                   ymin - ygap,ymax + ygap, prev_z + thickness,
                   areaBound, topverts, topfaces, top_orig);
}
