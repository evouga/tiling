#include "tinyxml2.h"
#include "Slice.h"
#include <sstream>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>

using namespace std;
using namespace tinyxml2;

void removeDuplicateVerts(Contour &c)
{
	vector<double> newx;
	vector<double> newy;
	int numpts = c.x.size();
	for(int i=0; i<numpts; i++)
	{
		int prev = (i==0 ? numpts-1 : i-1);
		if(c.x[prev] != c.x[i] && c.y[prev] != c.y[i])
		{
			newx.push_back(c.x[i]);
			newy.push_back(c.y[i]);
		}
	}
	c.x = newx;
	c.y = newy;
}

class ContourVisitor : public XMLVisitor
{
	public:
	ContourVisitor(const char *cname) : contourname(cname) {}

	virtual bool VisitEnter (const XMLElement &e, const XMLAttribute *a)
	{
		if(strcmp(e.Name(), "Contour"))
			return true;
		const char *cname = e.Attribute("name");
		if(!strcmp(cname, contourname))
		{
			const char *val = e.Attribute("points");
			stringstream vals(val);
			Contour contour;
			while(vals)
			{
				double x,y;
				char dummy;
				vals >> x;
				vals >> y;
				vals >> dummy;
				contour.x.push_back(x);
				contour.y.push_back(y);
			}
			removeDuplicateVerts(contour);
			contours.push_back(contour);
		}
		return true;
	}

	vector<Contour> contours;

	private:
	const char *contourname;
};

Slice *readSlice(const char *filename, const char *objectname)
{
	ContourVisitor cv(objectname);
	XMLDocument doc;

	if(doc.LoadFile(filename))
		return NULL;

	doc.Accept(&cv);

  double thickness;

	if(doc.FirstChildElement("Section")->QueryDoubleAttribute("thickness", &thickness) != XML_NO_ERROR) {
    printf("Returning NULL from %s:%d\n", __FILE__, __LINE__);
		return NULL;
	}

	return new Slice(cv.contours, thickness);
}

void readSlicesFromFolder(const char *baseFilename, const char *objectname, vector<Slice *> &slices) {
	slices.clear();
	int curslice = 0;
	while (true) {
		stringstream ss;
		ss << baseFilename;
		ss << ".";
		ss << curslice;
		Slice *slice = NULL;
		if ((slice = readSlice(ss.str().c_str(), objectname)))
			slices.push_back(slice);
		else
      break;
		curslice++;
	}
  // Give all contours unique ids.
  int contour_id = 10;
  for (Slice *slice : slices)
    for (Contour &contour : slice->contours)
      contour.contour_id = contour_id++;
  cout << "Slices: " << slices.size() << endl;
  cout << "Contours: " << contour_id << endl;
}

