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

Tile::Tile(const Slice &bottom, const Slice &top, bool scale) : bottom_(bottom), top_(top), use_scaling_(scale)
{
	computeCubeTransformation();
}

Tile::~Tile()
{
}

void Tile::computeCubeTransformation()
{
	Vector2d mins;
	Vector2d maxs;
	for(int i=0; i<2; i++)
	{
		mins[i] = numeric_limits<double>::infinity();
		maxs[i] = -numeric_limits<double>::infinity();
	}
	int numbottomcontours = bottom_.contours.size();
	for(int i=0; i<numbottomcontours; i++)
	{
		int numpts = bottom_.contours[i].x.size();
		for(int j=0; j<numpts; j++)
		{
			mins[0] = min(mins[0], bottom_.contours[i].x[j]);
			maxs[0] = max(maxs[0], bottom_.contours[i].x[j]);
			mins[1] = min(mins[1], bottom_.contours[i].y[j]);
			maxs[1] = max(maxs[1], bottom_.contours[i].y[j]);
		}
	}
	int numtopcontours = top_.contours.size();
	for(int i=0; i<numtopcontours; i++)
	{
		int numpts = top_.contours[i].x.size();
		for(int j=0; j<numpts; j++)
		{
			mins[0] = min(mins[0], top_.contours[i].x[j]);
			maxs[0] = max(maxs[0], top_.contours[i].x[j]);
			mins[1] = min(mins[1], top_.contours[i].y[j]);
			maxs[1] = max(maxs[1], top_.contours[i].y[j]);
		}
	}
  if (use_scaling_) {
    translate_ = -GLOBAL::z_lim*(mins+maxs);
    Vector2d widths =  (1.0 + tilePadding_)*(maxs-mins);
    scale_.setZero();
    scale_(0,0) = 1.0/widths[0];
    scale_(1,1) = 1.0/widths[1];
  } else {
    translate_.setZero();
    scale_.setZero();
    scale_(0,0) = 1;
    scale_(1,1) = 1;
  }
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

void Tile::scale(Vector2d &pt) {
	pt = scale_*(pt+translate_);
}

void Tile::unscale(Vector2d &pt) {
	pt = scale_.inverse()*pt - translate_;
}

// Adds original points to V and E. Should be already the correct size
int Tile::addOrig(const Slice &s, 
									Eigen::MatrixXd &V, Eigen::MatrixXi &E,
									Eigen::VectorXi &VM, Eigen::VectorXi &EM,
									int offset) {
	int numcontours = s.contours.size();
	for(int i=0; i<numcontours; i++)
	{
		int numpts = s.contours[i].x.size();
		for(int j=0; j<numpts; j++)
		{
			Vector2d pt(s.contours[i].x[j], s.contours[i].y[j]);
			//V.row(offset+j) = scale_*(pt+translate_);
			scale(pt);
			V(offset+j, 0) = pt(0);
			V(offset+j, 1) = pt(1);
			int prev = (j == 0 ? numpts-1 : j-1);
			E(offset+j, 0) = offset + j;
			E(offset+j, 1) = offset + prev;
			VM(offset+j) = GLOBAL::original_marker;
			EM(offset+j) = GLOBAL::original_marker;
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
	VM(0) = 0;
	VM(1) = 0;
	VM(2) = 0;
	VM(3) = 0;
	EM(0) = 0;
	EM(1) = 0;
	EM(2) = 0;
	EM(3) = 0;

	addOrig(s, V,E,VM,EM, 4);

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
  double xmin,xmax, ymin,ymax;
  if (use_scaling_) {
    xmin = ymin = -0.5;
    xmax = ymax = 0.5;
  } else {
    bottom_.xminmax(xmin,xmax);
    bottom_.yminmax(ymin,ymax);
    double txmin,txmax, tymin,tymax;
    top_.xminmax(txmin,txmax);
    top_.yminmax(tymin,tymax);
    xmin = min(xmin, txmin);
    xmax = max(xmax, txmax);
    ymin = min(ymin, tymin);
    ymax = max(ymax, tymax);
  }

  double thickness = bottom_.thickness;
#ifdef ZSCALE_HACK
  thickness = (std::max(ymax, xmax) - std::min(ymin, xmin)) / 2.0;
  printf("[%s:%d] Hack in place; thickness is set to %lf, instead of %lf\n",
         __FILE__, __LINE__, thickness, bottom_.thickness);
#endif
	triangulateSlice(bottom_, xmin,xmax, ymin,ymax, -thickness/2.0,
                   areaBound, botverts, botfaces, bot_orig);
	triangulateSlice(top_, xmin,xmax, ymin,ymax, thickness/2.0,
                   areaBound, topverts, topfaces, top_orig);
}
