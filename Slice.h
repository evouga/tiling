#ifndef SLICE_H
#define SLICE_H

#include <cstdio>
#include <vector>

#include "inpoly.hpp"

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

	bool isInside(double x, double y) const {
		for (unsigned int i = 0; i < contours.size(); ++i) {
			//printf("Checking for %lf,%lf inside index %d with size %lu:",
			//			 x, y, i, contours[i].x.size());
			std::vector<std::vector<double> > points = {contours[i].x, contours[i].y};
			if (inpoly(points, x, y)) {
				//printf(" true!\n");
				return true;
			} else {
				//printf(" false.\n");
			}
		}
		return false;
	}
};

#endif
