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
{ }

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
														Eigen::VectorXi &orig) {
	cout << "num contours: " << s.contours.size() << endl;

	int totpts = s.getNumPts();

  // Need alter vertices to include bounds.
  Eigen::MatrixXd V(totpts + 4, 2);
  Eigen::MatrixXi E(totpts + 4, 2);
  Eigen::MatrixXd H(0, 2);
  Eigen::VectorXi VM(totpts + 4);
  Eigen::VectorXi EM(totpts + 4);

  Eigen::MatrixXd lims;
	addOrig(s, V, E, VM, EM, lims, 0);

  // Add the bounding points after adding the original ones.
  V.row(totpts + 0) << xmin, ymin;
  V.row(totpts + 1) << xmax, ymin;
  V.row(totpts + 2) << xmax, ymax;
  V.row(totpts + 3) << xmin, ymax;

  // Attach new edges, label new edges and vertices.
  for (int i = 0; i < 4; ++i) {
    E(totpts + i, 0) = totpts + i;
    E(totpts + i, 1) = totpts + (i + 1) % 4;

    VM(totpts + i) = GLOBAL::nonoriginal_marker;
    EM(totpts + i) = GLOBAL::nonoriginal_marker;
  }

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

	cout << "Delaunay triangulating mesh with " << V.rows() << " verts" << endl;

	stringstream ss;
	ss << "Da" << areaBound << "q";
	MatrixXd V2;
	igl::triangle::triangulate(V, E, H, VM, EM, ss.str().c_str(), V2, faces, orig);
	cout << "done" << endl;
	verts.resize(V2.rows(), 3);
	flood_fill(faces, orig);

	for(int i=0; i<verts.rows(); i++) {
		verts(i, 0) = V2(i,0);
		verts(i, 1) = V2(i,1);
		verts(i, 2) = z;
	}
}

void Tile::triangulateSlices(double areaBound,
                             Eigen::MatrixXd &botV, Eigen::MatrixXi &botF,
                             Eigen::MatrixXd &topV, Eigen::MatrixXi &topF,
                             Eigen::VectorXi &botO, Eigen::VectorXi &topO) {
#ifdef ZSCALE_HACK
  double thickness = 0.25;
  // thickness = (botverts.colwise().maxCoeff() - botverts.colwise().minCoeff()).maxCoeff();

  printf("[%s:%d] Hack in place; thickness is set to %lf, instead of %lf\n",
         __FILE__, __LINE__, thickness, bottom_.thickness);
#else
  double thickness = bottom_.thickness;
#endif

  double xmin = bbox_(0, 0), xmax = bbox_(0, 1),
         ymin = bbox_(1, 0), ymax = bbox_(1, 1),
         zmin = bbox_(2, 0), zmax = zmin + thickness;

  cout << "Bounding box:" << endl;
  cout << "x: " << xmin << " " << xmax << endl;
  cout << "y: " << ymin << " " << ymax << endl;
  cout << "z: " << zmin << " " << zmax << endl;

  double xgap = (xmax - xmin) * GLOBAL::padding_perc,
         ygap = (ymax - ymin) * GLOBAL::padding_perc;


  double xMinPadded = xmin - xgap,
         xMaxPadded = xmax + xgap,
         yMinPadded = ymin - ygap,
         yMaxPadded = ymax + ygap;

  // If the bottom vertices is not empty, just use what you have.
  if (botV.rows() == 0) {
    triangulateSlice(bottom_,
                     xMinPadded, xMaxPadded,
                     yMinPadded, yMaxPadded, 0.0,
                     areaBound, botV, botF, botO);
  }

  double prev_z = botV(0, 2);
  triangulateSlice(top_,
                   xMinPadded, xMaxPadded,
                   yMinPadded, yMaxPadded, prev_z + thickness,
                   areaBound, topV, topF, topO);
}
