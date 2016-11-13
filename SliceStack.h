#ifndef SLICESTACK_H
#define SLICESTACK_H

#include <vector>
#include <iostream>

#include <Eigen/Core>

class Slice;

class SliceStack {
public:
	SliceStack(const char *baseFilename, const char *objectname);
	~SliceStack();

	int getNumSlices() {
    return numSlices_;
  }

	void triangulateSlice(int bottomidx, double areaBound,
                        Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces,
                        Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
                        Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig);

  // Triangulate the side, keeping a certain constant coordinate (see .cpp file)
  // All of these points will be along fixedCoord.
  void triangulateSide(int constantCoord,
                       double fixedCoord,
                       std::vector<Eigen::Vector3d> &verts,
                       Eigen::MatrixXd &V, Eigen::MatrixXi &F);

  // Will tetrahedralize slice and
  // return it as TV (vertices) TT (tets) and TF (faces)
  void tetrahedralizeSlice(const Eigen::MatrixXd &botverts,
                           const Eigen::MatrixXi &botfaces,
                           const Eigen::MatrixXd &topverts,
                           const Eigen::MatrixXi &topfaces,
                           const Eigen::VectorXi &botorig,
                           const Eigen::VectorXi &toporig,
                           Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                           Eigen::MatrixXi &TF, Eigen::VectorXi &TO);

	// Solve the laplacian and color by boundary conditions
	// TO is a vector corresponding to the vertices, where 1 means it's original
	// and 0 means it's not.
  // The output is Z, a vector that contains heat values for each vertex of TV.
	void computeLaplace(int slice_no, const Eigen::MatrixXd &TV,
											const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF,
											const Eigen::VectorXi &TO, Eigen::VectorXd &Z);

  int getSizeAt(int i);

private:
	void flipNormal(Eigen::MatrixXi &f);
	int numSlices_;
	std::vector<Slice*> slices_;
	std::vector<double> heights_;
  Eigen::MatrixXd bbox_;

  // Variable used by triangle to determine the maximum area of each triangle.
  // See -a flag. Lower value creates smaller triangles. With "hack" in place
  // (see Tile.cpp), a value of 0.04 isn't terrible.
  //const float triangle_max_area = 0.04;
  const float triangle_max_area = 0.01;
  // Variable used by tetgen to determine the maximum allowed *ratio* between
  // the radius and the area. See -q flag. With "hack" in place (see Tile.cpp),
  // a value of 1.4 isn't terrible.
  //const float tetgen_max_rad_ratio = 1.4;
  const float tetgen_max_rad_ratio = 1.1;
};

#endif
