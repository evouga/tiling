#ifndef SLICE_H
#define SLICE_H

#include <vector>

struct Contour
{
	std::vector<double> x;
	std::vector<double> y;
};

struct Slice
{
	std::vector<Contour> contours;
	double thickness;
};

#endif
