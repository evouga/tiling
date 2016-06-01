#ifndef SLICESTACK_H
#define SLICESTACK_H

#include <vector>
#include <iostream>

#include <Eigen/Core>

class Slice;

class SliceStack
{
public:
	SliceStack(const char *baseFilename, const char *objectname);
	~SliceStack();

	int getNumSlices() {return numSlices_;}

	void triangulateSlice(int bottomidx, double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
		Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig);

  void triangulateSide (int constantCoord, std::vector<Eigen::Vector3d> &verts,
    Eigen::MatrixXd &V, Eigen::MatrixXi &F);

	// Will tetrahedralize slice and return it as TV (vertices) TT (tets) and
	// TF (faces)
	void tetrahedralizeSlice(
			const Eigen::MatrixXd &botverts, const Eigen::MatrixXi &botfaces,
			const Eigen::MatrixXd &topverts, const Eigen::MatrixXi &topfaces,
			const Eigen::VectorXi &botorig, const Eigen::VectorXi &toporig,
			Eigen::MatrixXd &TV, Eigen::MatrixXi &TT, Eigen::MatrixXi &TF, Eigen::VectorXi &TO);

	// Solve the laplacian and color by boundary conditions
	// TO is a vector corresponding to the vertices, where 1 means it's original
	// and 0 means it's not.
	void computeLaplace(int slice_no, const Eigen::MatrixXd &TV,
											const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF,
											const Eigen::VectorXi &TO);

  int getSizeAt(int i);

private:
	void flipNormal(Eigen::MatrixXi &f);
	int numSlices_;
	std::vector<Slice *> slices_;
	std::vector<double> heights_;
};

#endif
