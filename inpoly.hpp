#ifndef INPOLY_HPP
#define INPOLY_HPP
#include <vector>

/**
 * Code adapted from http://www.visibone.com/inpoly/
 */
template<typename T>
bool                               /*  true=inside, false=outside          */
inpoly(                            /* is target point inside a 2D polygon? */
		std::vector<std::vector<T> > poly,            /*   polygon points, [0]=x, [1]=y       */
		T xt,                   /*   x (horizontal) of target point     */
		T yt,                   /*   y (vertical) of target point       */
		bool c_major = true)    /*   whether poly is stored in column-major order */
{
	int npoints;
	if (c_major)
		npoints = poly[0].size();
	else 
		npoints = poly.size();
	T xnew,ynew;
	T xold,yold;
	T x1,y1;
	T x2,y2;
	int i;
	bool inside = false;

	if (npoints < 3) {
		return(false);
	}
	if (c_major) {
		xold=poly[0][npoints-1];
		yold=poly[1][npoints-1];
	} else {
		xold=poly[npoints-1][0];
		yold=poly[npoints-1][1];
	}
	for (i=0 ; i < npoints ; i++) {
		if (c_major) {
			xnew=poly[0][i];
			ynew=poly[1][i];
		} else {
			xnew=poly[i][0];
			ynew=poly[i][1];
		}
		// Border points are missed
		if (xnew == xt && ynew == yt) {
			return true;
		}
		if (xnew > xold) {
			x1=xold;
			x2=xnew;
			y1=yold;
			y2=ynew;
		}
		else {
			x1=xnew;
			x2=xold;
			y1=ynew;
			y2=yold;
		}
		if ( ((xnew < xt) == (xt <= xold))          /* edge "open" at one end */
				&& ((yt-y1)*(x2-x1)
				    < (y2-y1)*(xt-x1)) ) {
			inside = !inside;
			printf("Switching!\n");
		}
		xold=xnew;
		yold=ynew;
	}
	return inside;
}
#endif // INPOLY_HPP
