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

	// for debugging purposes
	void triangulateSlice(int bottomidx, double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces);

	void tetrahedralizeSlice(const Eigen::MatrixXd &botverts,
    const Eigen::MatrixXi &botfaces, const Eigen::MatrixXd &topverts,
    const Eigen::MatrixXi &topfaces);

  int getSizeAt(int i);

private:
	int numSlices_;
	std::vector<Slice *> slices_;
	std::vector<double> heights_;
};

#endif
