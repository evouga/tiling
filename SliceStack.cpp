#include "SliceStack.h"
#include "SliceParser.h"
#include "Slice.h"
#include "Tile.h"

SliceStack::SliceStack(const char *baseFilename, const char *objectname)
{
	readSlicesFromFolder(baseFilename, objectname, slices_);
	numSlices_ = slices_.size();
	double z = 0;
	for(int i=0; i<numSlices_; i++)
	{
		heights_.push_back(z);
		z += slices_[i]->thickness;
	}
}

SliceStack::~SliceStack()
{
	for(int i=0; i<numSlices_; i++)
		delete slices_[i];
}

void SliceStack::triangulateSlice(int bottomidx, double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces)
{
	assert(bottomidx >= 0 && bottomidx < slices_.size());
	Tile t(*slices_[bottomidx], *slices_[bottomidx+1]);
	t.triangulateSlices(areaBound, botverts, botfaces, topverts, topfaces);
}

