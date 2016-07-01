#ifndef SLICE_H
#define SLICE_H

#include <algorithm>
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

  void xminmax(double &mn, double &mx) const {
    bool first = true;
    for (unsigned int i = 0; i < contours.size(); ++i) {
      auto el = std::minmax_element(contours[i].x.begin(), contours[i].x.end());
      if (first) {
        first = false;
        mn = *el.first;
        mx = *el.second;
      } else {
        mn = std::min(mn, *el.first);
        mx = std::max(mx, *el.second);
      }
    }
  }
  void yminmax(double &mn, double &mx) const {
    bool first = true;
    for (unsigned int i = 0; i < contours.size(); ++i) {
      auto el = std::minmax_element(contours[i].y.begin(), contours[i].y.end());
      if (first) {
        first = false;
        mn = *el.first;
        mx = *el.second;
      } else {
        mn = std::min(mn, *el.first);
        mx = std::max(mx, *el.second);
      }
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
