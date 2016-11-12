#ifndef SLICE_H
#define SLICE_H

#include <algorithm>
#include <cstdio>
#include <limits>
#include <vector>

//#include "inpoly.hpp"

struct Contour {
	std::vector<double> x;
	std::vector<double> y;
};

struct Slice {
  std::vector<Contour> contours;
	double thickness;
  double minX, minY, maxX, maxY;

  Slice(const std::vector<Contour>& contours, double thickness) {
    this->contours = contours;
    this->thickness = thickness;

    this->minX = std::numeric_limits<double>::max();
    this->maxX = std::numeric_limits<double>::min();
    this->minY = std::numeric_limits<double>::max();
    this->maxY = std::numeric_limits<double>::min();

    for (unsigned int i = 0; i < contours.size(); ++i) {
      auto elX = std::minmax_element(contours[i].x.begin(), contours[i].x.end());
      this->minX = std::min(minX, *elX.first);
      this->maxX = std::max(maxX, *elX.second);

      auto elY = std::minmax_element(contours[i].y.begin(), contours[i].y.end());
      this->minY = std::min(minY, *elY.first);
      this->maxY = std::max(maxY, *elY.second);
    }
  }

	unsigned int getNumPts() const {
		unsigned int totpts = 0;
		for (unsigned int i = 0; i < contours.size(); ++i) {
			totpts += contours[i].x.size();
		}
		return totpts;
	}
};

#endif
