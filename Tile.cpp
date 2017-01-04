#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include <igl/jet.h>
#include <igl/triangle/triangulate.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "Tile.h"
#include "Slice.h"

using namespace Eigen;
using namespace std;

#define ZSCALE_HACK

const double Tile::tilePadding_ = 0.10;

Tile::Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox)
  : bottom_(bottom), top_(top), bbox_(bbox) {
    vector<int> empty;
  allowed_bot_ = empty;
  allowed_top_ = empty;
}

Tile::Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox,
           const std::vector<int> &allowed_bot, const std::vector<int> &allowed_top)
: bottom_(bottom), top_(top), bbox_(bbox),
    allowed_bot_(allowed_bot), allowed_top_(allowed_top)
{ }

Tile::~Tile()
{ }

// Adds original points to V and E. Should be already the correct size
void Tile::addOrig(const Slice &s,
                   std::vector<RowVector2d> &V,
                   std::vector<RowVector2i> &E,
                   std::vector<int> &VM, std::vector<int> &EM,
                   MatrixXd &lims,
                   const std::vector<int> &allowed) {
  bool first = true;
  lims.resize(2, 2);
  int numcontours = s.contours.size();
  // Iterate over all contours
  for(int i=0; i<numcontours; i++) {
    // Check and make sure this is a valid contour.
    bool is_allowed = false;
    for (int j = 0; j < allowed.size(); ++j) {
      if (i == allowed[j]) {
        is_allowed = true;
        break;
      }
    }
    // If allowed is empty, all contours are valid.
    if (!is_allowed && allowed.size() > 0) {
      continue;
    }

    // Get the number of points for this contour.
    int numpts = s.contours[i].x.size();

    // Number of points we've added to this point.
    int offset = V.size();

    // Find all points in this contour.
    for(int j=0; j<numpts; j++) {
      // pt is the point.
      RowVector2d pt(s.contours[i].x[j], s.contours[i].y[j]);
      V.push_back(pt);

      // Which point should we connect this to?
      int prev = (j == 0 ? numpts-1 : j-1);
      // ed is the edge.
      RowVector2i ed(offset + j, offset + prev);
      E.push_back(ed);

      // The vertex will be marked by the contour id.
      VM.push_back(s.contours[i].contour_id);
      EM.push_back(GLOBAL::original_marker);

      if (first) {
        lims(0, 0) = lims(0, 1) = pt(0);
        lims(1, 0) = lims(1, 1) = pt(1);
        first = false;
      } else {
        lims(0, 0) = min(lims(0, 0), pt(0));
        lims(0, 1) = max(lims(0, 1), pt(0));
        lims(1, 0) = min(lims(1, 0), pt(1));
        lims(1, 1) = max(lims(1, 1), pt(1));
      }
    }
  }
}

void flood_fill(const Eigen::MatrixXi &faces, Eigen::VectorXi &orig) {
  bool dirty = true;
  while(dirty) {
    dirty = false;
    for (int i = 0; i < faces.rows(); ++i) {
      if (orig(faces(i,0)) == 1 || orig(faces(i,1)) == 1 || orig(faces(i,2)) == 1) {
        for (int j = 0; j < 3; ++j) {
          int vert = faces(i,j);
          if (orig(vert) == 0) {
            orig(vert) = 1;
            dirty = true;
          }
        }
      }
    }
  }

  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) == 0) {
      orig(i) = GLOBAL::original_marker;
    }
  }

  dirty = true;
  while(dirty) {
    dirty = false;
    for (int i = 0; i < faces.rows(); ++i) {
      int marker = GLOBAL::nonoriginal_marker;
      for (int j = 0; j < 3; j++)
        marker = max(marker, orig(faces(i, j)));
      for (int j = 0; j < 3; ++j) {
        int vert = faces(i,j);
        if (orig(vert) >= GLOBAL::original_marker && orig(vert) != marker) {
          orig(vert) = marker;
          dirty = true;
        }
      }
    }
  }
}

// Triangulates a single contour.
void Tile::triangulateSlice(const Slice &s, double z, double areaBound,
                            Eigen::MatrixXd &verts, Eigen::MatrixXi &faces,
                            Eigen::VectorXi &orig,
                            const std::vector<int> &allowed) {
  cout << "num contours: " << s.contours.size() << endl;

  std::vector<RowVector2d> V_v;
  std::vector<RowVector2i> E_v;
  std::vector<int> VM_v, EM_v;
  Eigen::MatrixXd lims;
  addOrig(s, V_v, E_v, VM_v, EM_v, lims, allowed);

  // Need alter vertices to include bounds.
  int totpts = V_v.size();
  Eigen::MatrixXd V(totpts + 4, 2);
  Eigen::MatrixXi E(totpts + 4, 2);
  Eigen::MatrixXd H(0, 2);
  Eigen::VectorXi VM(totpts + 4);
  Eigen::VectorXi EM(totpts + 4);

  // Copy stuff from vectors. They all have the same size.
  for (int i = 0; i < totpts; ++i) {
    V.row(i) = V_v[i];
    E.row(i) = E_v[i];
    VM(i) = VM_v[i];
    EM(i) = EM_v[i];
  }

  // Add the bounding points after adding the original ones.
  double xgap = (s.maxX - s.minX) * GLOBAL::padding_perc,
         ygap = (s.maxY - s.minY) * GLOBAL::padding_perc;
  V.row(totpts + 0) << s.minX - xgap, s.minY - ygap;
  V.row(totpts + 1) << s.maxX + xgap, s.minY - ygap;
  V.row(totpts + 2) << s.maxX + xgap, s.maxY + ygap;
  V.row(totpts + 3) << s.minX - xgap, s.maxY + ygap;

  // Attach new edges, label new edges and vertices.
  for (int i = 0; i < 4; ++i) {
    E(totpts + i, 0) = totpts + i;
    E(totpts + i, 1) = totpts + (i + 1) % 4;

    VM(totpts + i) = GLOBAL::nonoriginal_marker;
    EM(totpts + i) = GLOBAL::nonoriginal_marker;
  }

  cout << "Delaunay triangulating mesh with " << V.rows() << " verts" << endl;

  stringstream ss;
  ss << "QDa" << areaBound << "q";

  MatrixXd V2;
  igl::triangle::triangulate(V, E, H, VM, EM, ss.str().c_str(), V2, faces, orig);

  verts.resize(V2.rows(), 3);
  flood_fill(faces, orig);

  // 3D-fying the tile.
  for(int i=0; i<verts.rows(); i++) {
    verts(i, 0) = V2(i, 0);
    verts(i, 1) = V2(i, 1);
    verts(i, 2) = z;
  }
}

void Tile::triangulateSlices(double areaBound,
                             Eigen::MatrixXd &botV, Eigen::MatrixXi &botF,
                             Eigen::MatrixXd &topV, Eigen::MatrixXi &topF,
                             Eigen::VectorXi &botO, Eigen::VectorXi &topO) {
#ifdef ZSCALE_HACK
  double thickness = 0.4;
  // thickness = (botverts.colwise().maxCoeff() - botverts.colwise().minCoeff()).maxCoeff();
  printf("[%s:%d] Hack in place; thickness is set to %lf, instead of %lf\n",
         __FILE__, __LINE__, thickness, bottom_.thickness);
#else
  double thickness = bottom_.thickness;
#endif

  // If the bottom vertices is not empty, just use what you have.
  if (botV.rows() == 0) {
    triangulateSlice(bottom_, 0.0, areaBound, botV, botF, botO, allowed_bot_);
  }

  triangulateSlice(top_, botV(0, 2) + thickness, areaBound, topV, topF, topO, allowed_top_);
}
