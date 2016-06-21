#ifndef SLICE_H
#define SLICE_H

#include <cstdio>
#include <vector>

//#include "inpoly.hpp"

struct Contour
{
	std::vector<double> x;
	std::vector<double> y;
};

struct Slice
{
	std::vector<Contour> contours;
	double thickness;
	unsigned int getNumPts() const {
		unsigned int totpts = 0;
		for (unsigned int i = 0; i < contours.size(); ++i) {
			totpts += contours[i].x.size();
		}
		return totpts;
	}
};

#endif
