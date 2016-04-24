#include "Tile.h"
#include "Slice.h"
#include <limits>
#include <igl/triangle/triangulate.h>
#include <sstream>
#include <iostream>

using namespace Eigen;
using namespace std;

const double Tile::tilePadding_ = 0.10;

Tile::Tile(const Slice &bottom, const Slice &top) : bottom_(bottom), top_(top)
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
	translate_ = -0.5*(mins+maxs);
	Vector2d widths =  (1.0 + tilePadding_)*(maxs-mins);
	scale_.setZero();
	scale_(0,0) = 1.0/widths[0];
	scale_(1,1) = 1.0/widths[1];
}	

void Tile::triangulateSlice(const Slice &s, double z, double areaBound, Eigen::MatrixXd &verts, Eigen::MatrixXi &faces)
{
	int numcontours = s.contours.size();
	cout << "num contours: " << numcontours << endl;
	int totpts = 0;
	for(int i=0; i<numcontours; i++)
	{
		int numpts = s.contours[i].x.size();
		totpts += numpts;
	}

	MatrixXd V(totpts+4, 2);
	MatrixXi E(totpts+4, 2);
	MatrixXd H(0, 2);

	V(0,0) = -1.0;
	V(0,1) = -1.0;
	V(1,0) = 1.0;
	V(1,1) = -1.0;
	V(2,0) = 1.0;
	V(2,1) = 1.0;
	V(3,0) = -1.0;
	V(3,1) = 1.0;
	E(0,0) = 0;
	E(0,1) = 1;
	E(1,0) = 1;
	E(1,1) = 2;
	E(2,0) = 2;
	E(2,1) = 3;
	E(3,0) = 3;
	E(3,1) = 0;

	int offset = 4;

	for(int i=0; i<numcontours; i++)
	{
		int numpts = s.contours[i].x.size();
		for(int j=0; j<numpts; j++)
		{
			Vector2d pt(s.contours[i].x[j], s.contours[i].y[j]);
			V.row(offset+j) = scale_*(pt+translate_);
			int prev = (j == 0 ? numpts-1 : j-1);
			E(offset+j, 0) = offset + j;
			E(offset+j, 1) = offset + prev;
		}
		offset += numpts;
	}
	stringstream ss;
	ss << "Da" << areaBound << "q";

	MatrixXd V2;
	cout << "Delaunay triangulating mesh with " << V.rows() << " verts" << endl;
	igl::triangle::triangulate(V,E,H,ss.str().c_str(),V2,faces);	
	cout << "done" << endl;
	verts.resize(V2.rows(), 3);
	for(int i=0; i<verts.rows(); i++)
	{
		verts(i, 0) = V2(i,0);
		verts(i, 1) = V2(i,1);
		verts(i, 2) = z;
	}
}

void Tile::triangulateSlices(double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces)
{
	triangulateSlice(bottom_, -1.0, areaBound, botverts, botfaces);
	triangulateSlice(top_, 1.0, areaBound, topverts, topfaces);
}
