// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Nathan Clement <nathanlclement@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <map>
#include <vector>

#include <Eigen/Core>

#include <igl/bfs_orient.h>
#include <igl/collapse_small_triangles.h>
#include <igl/remove_unreferenced.h>
#include <igl/orientable_patches.h>
#include <igl/orient_outward.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_duplicates.h>
#include <igl/resolve_duplicated_faces.h>

#include "glob_defs.h"
#include "Helpers.h"

namespace {
int oTimes[5];
struct MarchingTetsDat {
  const Eigen::MatrixXd &V;
  const Eigen::VectorXd &H;
  std::vector<Eigen::RowVector3i> &faces;
  std::vector<int> &faces_markers;
  const std::map<std::vector<int>, int> &face_counts;
  std::vector<Eigen::RowVector3d> &new_verts;
  double offset;

  MarchingTetsDat(const Eigen::MatrixXd &_V, const Eigen::VectorXd &_H,
                  std::vector<Eigen::RowVector3i> &_f,
                  std::vector<int> &_fm,
                  const std::map<std::vector<int>, int> &_fc,
                  std::vector<Eigen::RowVector3d> &_nv,
                  double _off)
      : V(_V), H(_H), faces(_f), faces_markers(_fm), face_counts(_fc),
        new_verts(_nv), offset(_off)
  {}

};

// Returns a new point between a and b. how_much is the distance between the
// two, where the result will be (how_much * a + (1 - how_much) * b) / 2
//
// Assume that aval is <= offset and bval >= offset.
Eigen::RowVector3d split(const Eigen::RowVectorXd &a,
                         const Eigen::RowVectorXd &b,
                         double aval, double bval,
                         double offset) {
  Eigen::RowVector3d newP;
  double how_much = (offset - aval) / (bval - aval);
  if (how_much < 0 || how_much > 1) {
    fprintf(stderr, "Error: how_much is %lf (from %lf %lf %lf)\n",
            how_much, aval, bval, offset);
  }
  return how_much * a + (1 - how_much) * b;
}


int getMapCounts(const std::map<std::vector<int>, int> &counts,
                 std::vector<int> &c) {
  std::sort(c.begin(), c.end());
  if (counts.find(c) == counts.end()) {
    printf("Error!! Couldn't find it!!\n");
  }
  return counts.find(c)->second;
}
int getMapCounts(const std::map<std::vector<int>, int> &counts,
                 const Eigen::RowVector3i &f) {
  std::vector<int> t({f(0), f(1), f(2)});
  return getMapCounts(counts, t);
}
bool isValid(std::vector<int> &c,
             const std::map<std::vector<int>, int> &counts,
             const Eigen::VectorXd &H) {
  return getMapCounts(counts, c) == 1;
}
bool isValid(const Eigen::RowVector3i &f,
             const std::map<std::vector<int>, int> &counts,
             const Eigen::VectorXd &H,
             bool test_temp = false) {
  if (test_temp) {
    for (int i = 0; i < 3; ++i) {
      if (H(f(i)) == GLOBAL::outside_temp) return false;
    }
  }

  return getMapCounts(counts, f) == 1;
}

void addOrig(const Eigen::RowVectorXi &T, int marker,
             std::vector<Eigen::RowVector3i> &faces,
             std::vector<int> &faces_markers) {
  Eigen::RowVector3i f1, f2, f3, f4;
  f1 << T(0), T(1), T(3);
  f2 << T(0), T(2), T(1);
  f3 << T(3), T(2), T(0);
  f4 << T(1), T(2), T(3);
  faces.push_back(f1);
  faces.push_back(f2);
  faces.push_back(f3);
  faces.push_back(f4);
  faces_markers.push_back(marker);
  faces_markers.push_back(marker);
  faces_markers.push_back(marker);
  faces_markers.push_back(marker);
}


// Will check the corner heat values, and return:
//  true if faces have been added; finished
//  false if this doesn't apply; add normally
bool heatCheck(MarchingTetsDat &dd,
               const std::vector<int> &inside, const std::vector<int> &outside,
               const std::vector<int> &inside_t,
               const std::vector<int> &outside_t,
               const Eigen::RowVectorXi &T,
               int marker) {
  // Turn off heat check
  return false;

  // A couple of cases to be concerned about with temperature.
  if (inside_t.size() > 0 && outside_t.size() > 0) {
    //return true;
    //addOrig(T, 5, dd.faces, dd.faces_markers); return true;

    // Impossible cases:
    //  inside_t: 1 outside_t: 3
    //  inside_t: 3 outside_t: 1
    //  inside_t: 2 outside_t: 2
    //
    // Possible cases:
    //  inside_t: 1 outside_t: 2 (4 faces)
    //  inside_t: 2 outside_t: 1 (3 faces)
    //  inside_t: 1 outside_t: 1 (5 faces)
    if (inside_t.size() == 1 && outside_t.size() == 2) {
      // If all three vertices are outside, this becomes just a line; ignore.
      if (outside.size() == 3) {
        return true;
      }

      // inside.size() == 2
      // Get the non-boundary vertex.
      int inV1 = T(inside_t[0]), inV2 = T(inside[0]);
      if (inV2 == inV1) inV2 = T(inside[1]);
      int outV1 = T(outside_t[0]), outV2 = T(outside_t[1]);

      // Create two new points.
      int newV_pos = dd.V.rows() + dd.new_verts.size();
      dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV1),
                                   dd.H(inV2), dd.H(outV1), dd.offset));
      dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV2),
                                   dd.H(inV2), dd.H(outV2), dd.offset));
      int h1 = newV_pos, h2 = newV_pos + 1;
      // If the outside face is on the boundary, add it.
      if (isValid({inV2, outV1, outV2}, dd.face_counts, dd.H)) {
        Eigen::RowVector3i newF({h1, h2, inV2});
        dd.faces.push_back(newF);
        dd.faces_markers.push_back(marker);
      }
      // Check for the right and left face as well.
      if (isValid({inV1, inV2, outV1}, dd.face_counts, dd.H)) {
        Eigen::RowVector3i newF({inV1, h1, inV2});
        dd.faces.push_back(newF);
        dd.faces_markers.push_back(marker);
      }
      if (isValid({inV1, inV2, outV2}, dd.face_counts, dd.H)) {
        Eigen::RowVector3i newF({inV1, h2, inV2});
        dd.faces.push_back(newF);
        dd.faces_markers.push_back(marker);
      }
      // Bottom, always add.
      Eigen::RowVector3i newF({inV1, h1, h2});
      dd.faces.push_back(newF);
      dd.faces_markers.push_back(marker);
      // Finished.
      return true;
    } else if (inside_t.size() == 2 && outside_t.size() == 1) {
      // Two cases:
      //  - inside:2 (both inside_t), outside: 2 (one outside_t)
      //  - inside:3 (two inside_t), outside: 1 (outside_t)
      if (inside.size() == 2) {
        // Split face into two faces if it's valid.
        // Create two new points.
        int inV1 = T(inside_t[0]), inV2 = T(inside_t[1]);
        int outV = outside[0] == outside_t[0] ? T(outside[1]) : T(outside[0]);

        // Want to always add these faces
        int newV_pos = dd.V.rows() + dd.new_verts.size();
        dd.new_verts.push_back(split(dd.V.row(inV1), dd.V.row(outV),
                                     dd.H(inV1), dd.H(outV), dd.offset));
        dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV),
                                     dd.H(inV2), dd.H(outV), dd.offset));
        int h1 = newV_pos, h2 = newV_pos + 1;

        Eigen::RowVector3i newF1, newF2;
        newF1 << inV1, h1, inV2;
        newF2 << inV2, h1, h2;
        dd.faces.push_back(newF1);
        dd.faces.push_back(newF2);
        dd.faces_markers.push_back(marker);
        dd.faces_markers.push_back(marker);
        return true;
      } else if (inside.size() == 3) {
        //  - inside:3 (two inside_t), outside: 1 (outside_t)
        // Three faces.
        int inV1 = T(inside_t[0]), inV2 = T(inside_t[1]);
        int inV3;
        for (int i = 0; i < inside.size(); ++i) {
          inV3 = T(inside[i]);
          if (inV3 != inV1 && inV3 != inV2) {
            break;
          }
        }
        // TODO here

        int outV = T(outside[0]);

        // Create one new point.
        int h = dd.V.rows() + dd.new_verts.size();
        dd.new_verts.push_back(split(dd.V.row(inV3), dd.V.row(outV),
                                     dd.H(inV3), dd.H(outV), dd.offset));
        // Left and right
        if (isValid({inV1, outV, inV3}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, h, inV3});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        if (isValid({inV2, outV, inV3}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV2, inV3, h});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Front
        if (isValid({inV1, inV2, inV3}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, inV2, inV3});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Bottom, always add.
        if (1) {
          printf("Diff is %lf vs %lf (%d vs %d)\n",
                 dd.H(inV3), dd.H(outV), outside[0], outside_t[0]);
          if (dd.V(inV1, 2) == dd.new_verts[h - dd.V.rows()][2]) {
            printf("Error: Parallel!!\n");
          }
          Eigen::RowVector3i newF({inV1, inV2, h});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        return true;
      } else {
        printf("[%s:%d] Warning: Unknown tet structure!!\n", __FILE__, __LINE__);
      }
      
    } else if (inside_t.size() == 1 && outside_t.size() == 1) {
      // Possible add 5 faces.
      // Two cases:
      //  - inside:3 (one inside_t), outside: 1 (outside_t): 5 faces
      //  - inside:2 (one inside_t), outside: 2 (one outside_t): 3 faces
      //  - inside:1 (inside_t), outside: 3 (one outside_t): 1 face
      if (inside.size() == 1) {
        int inV = T(inside[0]);
        int outV1, outV2;
        for (int i = 0; i < outside.size(); ++i) {
          if (outside[i] == outside_t[0]) {
            outV1 = T(outside[(i + 1) % 3]);
            outV2 = T(outside[(i + 2) % 3]);
            break;
          }
        }
        if (isValid({inV, outV1, outV2}, dd.face_counts, dd.H)) {
          // add two new verts and one new face.
          int newV_pos = dd.V.rows() + dd.new_verts.size();
          dd.new_verts.push_back(split(dd.V.row(inV), dd.V.row(outV1),
                                       dd.H(inV), dd.H(outV1), dd.offset));
          dd.new_verts.push_back(split(dd.V.row(inV), dd.V.row(outV2),
                                       dd.H(inV), dd.H(outV2), dd.offset));
          int h1 = newV_pos, h2 = newV_pos + 1;
          Eigen::RowVector3i newF({inV, h1, h2});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        return true;
      } else if (inside.size() == 2) {
        //  - inside:2 (one inside_t), outside: 2 (one outside_t): 3 faces
        // 2 new verts and 4 new faces.
        int inV1 = T(inside_t[0]);
        int inV2 = T(inside[0]) == inV1 ? T(inside[1]) : T(inside[0]);
        int outV1 = T(outside[0]), outV2 = T(outside[1]);
        int newV_pos = dd.V.rows() + dd.new_verts.size();
        dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV1),
                                     dd.H(inV2), dd.H(outV1), dd.offset));
        dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV2),
                                     dd.H(inV2), dd.H(outV2), dd.offset));
        int h1 = newV_pos, h2 = newV_pos + 1;
        // Add back face
        if (isValid({inV2, outV1, outV2}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV2, h1, h2});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Add side faces
        if (isValid({inV1, inV2, outV1}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, h1, inV2});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        if (isValid({inV1, inV2, outV2}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, inV2, h2});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Always add bottom face.
        Eigen::RowVector3i newF({inV1, h1, h2});
        dd.faces.push_back(newF);
        dd.faces_markers.push_back(marker);
        return true;
      } else if (inside.size() == 3) {
        //  - inside:3 (one inside_t), outside: 1 (outside_t): 5 faces
        // Create 2 new verts, 5 new faces
        int inV1 = T(inside_t[0]), inV2, inV3;
        for (int i = 0; i < inside.size(); ++i) {
          if (T(inside[i]) == inV1) {
            inV2 = T(inside[(i + 1) % 3]);
            inV3 = T(inside[(i + 2) % 3]);
            break;
          }
        }
        int outV = T(outside[0]);
        int newV_pos = dd.V.rows() + dd.new_verts.size();
        dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV),
                                     dd.H(inV2), dd.H(outV), dd.offset));
        dd.new_verts.push_back(split(dd.V.row(inV3), dd.V.row(outV),
                                     dd.H(inV3), dd.H(outV), dd.offset));
        int h2 = newV_pos, h3 = newV_pos + 1;

        // Front face
        if (isValid({inV2, inV3, outV}, dd.face_counts, dd.H)) {
          // two faces
          Eigen::RowVector3i newF1({inV2, h2, inV3});
          Eigen::RowVector3i newF2({h2, h3, inV3});
          dd.faces.push_back(newF1);
          dd.faces.push_back(newF2);
          dd.faces_markers.push_back(marker);
          dd.faces_markers.push_back(marker);
        }
        // Side faces.
        if (isValid({inV1, inV2, outV}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, h2, inV2});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        if (isValid({inV1, inV3, outV}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, h3, inV3});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Back (inside) face
        if (isValid({inV1, inV2, inV3}, dd.face_counts, dd.H)) {
          Eigen::RowVector3i newF({inV1, inV2, inV3});
          dd.faces.push_back(newF);
          dd.faces_markers.push_back(marker);
        }
        // Bottom face (always add)
        Eigen::RowVector3i newF({inV1, h2, h3});
        dd.faces.push_back(newF);
        dd.faces_markers.push_back(marker);
        return true;
      } else {
        printf("[%s:%d] Warning: Unknown tet structure!!\n", __FILE__, __LINE__);
      }
    } else {
      printf("[%s:%d] Warning: Unknown tet structure!!\n", __FILE__, __LINE__);
    }
  }

  // Couldn't find something useful.
  return false;
}

// (A) Takes care of:
//   inside: 1, identical: 3 (remove three duplicated faces later)
//   inside: 2, identical: 2 (remove all four duplicated faces later)
//   inside: 3, identical: 1 (remove all four duplicated faces later)
//   inside: 4, identical: 0 (remove all four duplicated faces later)
void case0Out(MarchingTetsDat &dd, const Eigen::RowVectorXi &T,
              const std::vector<int> &inside, const std::vector<int> &outside,
              const std::vector<int> &identical,
              const std::vector<int> &inside_t, const std::vector<int> &outside_t) {
  // Check the heat thing.
  if(heatCheck(dd, inside,outside, inside_t,outside_t, T, 0)) return;

  // Add all 4 faces.
  Eigen::RowVector3i f1, f2, f3, f4;
  f1 << T(0), T(1), T(3);
  f2 << T(0), T(2), T(1);
  f3 << T(3), T(2), T(0);
  f4 << T(1), T(2), T(3);
  if (isValid(f1, dd.face_counts, dd.H, true)) {
    dd.faces.push_back(f1);
    dd.faces_markers.push_back(0);
  }
  if (isValid(f2, dd.face_counts, dd.H, true)) {
    dd.faces.push_back(f2);
    dd.faces_markers.push_back(0);
  }
  if (isValid(f3, dd.face_counts, dd.H, true)) {
    dd.faces.push_back(f3);
    dd.faces_markers.push_back(0);
  }
  if (isValid(f4, dd.face_counts, dd.H, true)) {
    dd.faces.push_back(f4);
    dd.faces_markers.push_back(0);
  }
}

// (B) Takes care of:
//   inside: 1 outside: 3
//   inside: 1 outside: 2 identical: 1
//   inside: 1 outside: 1 identical: 2
void case1In(MarchingTetsDat &dd, const Eigen::RowVectorXi &T,
             const std::vector<int> &inside, const std::vector<int> &outside,
             const std::vector<int> &identical,
             const std::vector<int> &inside_t, const std::vector<int> &outside_t) {
  
  // Check the heat thing.
  if(heatCheck(dd, inside,outside, inside_t,outside_t, T, 2)) return;

  // Only one inside, need to divide each of the remaining three in half,
  // then add one new face (ignore the other three)
  int inV = T(inside[0]);
  int outV1, outV2, outV3;
  if (outside.size() == 3) {
    outV1 = T(outside[0]);
    outV2 = T(outside[1]);
    outV3 = T(outside[2]);
  } else if (outside.size() == 2) {
    outV1 = T(outside[0]);
    outV2 = T(outside[1]);
    outV3 = T(identical[0]);
  } else if (outside.size() == 1) {
    outV1 = T(outside[0]);
    outV2 = T(identical[0]);
    outV3 = T(identical[1]);
  }

  int newV_pos = dd.V.rows() + dd.new_verts.size();
  dd.new_verts.push_back(split(dd.V.row(inV), dd.V.row(outV1),
                            dd.H(inV), dd.H(outV1), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV), dd.V.row(outV2),
                            dd.H(inV), dd.H(outV2), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV), dd.V.row(outV3),
                            dd.H(inV), dd.H(outV3), dd.offset));

  // Create the new faces.
  Eigen::RowVector3i newF1, newF2, newF3, newF4;
  // Don't care about normals here... Can fix later (?)
  newF1 << newV_pos, newV_pos + 1, newV_pos + 2;
  newF2 << inV, newV_pos, newV_pos + 1;
  newF3 << inV, newV_pos + 1, newV_pos + 2;
  newF4 << inV, newV_pos + 2, newV_pos;
  // This face is automatically added.
  dd.faces.push_back(newF1);
  dd.faces_markers.push_back(2);
  // Only add the rest of the faces if they're on a a boundary.
  if (isValid({inV, outV1, outV2}, dd.face_counts, dd.H)) {
    dd.faces.push_back(newF2);
    dd.faces_markers.push_back(2);
  }
  if (isValid({inV, outV2, outV3}, dd.face_counts, dd.H)) {
    dd.faces.push_back(newF3);
    dd.faces_markers.push_back(2);
  }
  if (isValid({inV, outV3, outV1}, dd.face_counts, dd.H)) {
    dd.faces.push_back(newF4);
    dd.faces_markers.push_back(2);
  }
}

// (D) takes care of:
//   inside: 2 outside: 1 identical: 1
//   inside: 2 outside: 2 identical: 0
void case2In2Out(MarchingTetsDat &dd, const Eigen::RowVectorXi &T,
                 const std::vector<int> &inside, const std::vector<int> &outside,
                 const std::vector<int> &identical,
                 const std::vector<int> &inside_t, const std::vector<int> &outside_t) {
  // 2 inside, 2 outside. Insert four faces between 4 points.

  // Check the heat thing.
  if(heatCheck(dd, inside,outside, inside_t,outside_t, T, 3)) return;

  int outV1 = T(outside[0]);
  int outV2 = T(outside[1]);
  int inV1 = T(inside[0]);
  int inV2 = T(inside[1]);

  // Create 4 new points.
  int newV_pos = dd.V.rows() + dd.new_verts.size();
  dd.new_verts.push_back(split(dd.V.row(inV1), dd.V.row(outV1), // h11
                            dd.H(inV1), dd.H(outV1), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV1), // h21
                            dd.H(inV2), dd.H(outV1), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV1), dd.V.row(outV2), // h12
                            dd.H(inV1), dd.H(outV2), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV2), // h22
                            dd.H(inV2), dd.H(outV2), dd.offset));

  // Create two faces
  std::vector<std::vector<int> > newFs;

  // Special case:
  // Front: i2, h22, h21
  if (isValid({inV2, outV2, outV1}, dd.face_counts, dd.H)) {
    newFs.push_back({inV2, newV_pos + 3, newV_pos + 1});
  }
  // Back: i1, h11, h12
  if (isValid({inV1, outV1, outV2}, dd.face_counts, dd.H)) {
    newFs.push_back({inV1, newV_pos, newV_pos + 2});
  }
  // Side: i2, i1, h22 + h12, i1, h22
  if (isValid({inV1, inV2, outV2}, dd.face_counts, dd.H)) {
    newFs.push_back({inV2,         inV1, newV_pos + 3});
    newFs.push_back({newV_pos + 3, inV1, newV_pos + 2});
  }
  // Side: i2, i1, h21 + h21, i1, h11
  if (isValid({inV1, inV2, outV1}, dd.face_counts, dd.H)) {
    newFs.push_back({inV2,         inV1, newV_pos + 1});
    newFs.push_back({newV_pos + 1, inV1, newV_pos});
  }
  // Always add these.
  // Check that none of these are on the boundary.
  //if (H(outV1) != GLOBAL::outside_temp && H(outV2) != GLOBAL::outside_temp) {
  // Top: h22, h21, h12 + h12, h21, h11
  newFs.push_back({newV_pos + 3, newV_pos + 1, newV_pos + 2});
  newFs.push_back({newV_pos + 2, newV_pos + 1, newV_pos});
  //}

  for (int i = 0; i < newFs.size(); ++i) {
    Eigen::RowVector3i t({newFs[i][0], newFs[i][1], newFs[i][2]});
    dd.faces.push_back(t);
    dd.faces_markers.push_back(4);
  }
}

// (C) takes care of:
//   inside: 3 outside: 1
void case3In1Out(MarchingTetsDat &dd, const Eigen::RowVectorXi &T,
                 const std::vector<int> &inside, const std::vector<int> &outside,
                 const std::vector<int> &identical,
                 const std::vector<int> &inside_t, const std::vector<int> &outside_t) {
  // Only one outside, need to divide each of the remaining three in half,
  // then add one new face.

  // Check the heat thing.
  if(heatCheck(dd, inside,outside, inside_t,outside_t, T, 3)) return;

  int outV = T(outside[0]);
  int inV1 = T(inside[0]);
  int inV2 = T(inside[1]);
  int inV3 = T(inside[2]);

  int newV_pos = dd.V.rows() + dd.new_verts.size();
  dd.new_verts.push_back(split(dd.V.row(inV1), dd.V.row(outV),
                            dd.H(inV1), dd.H(outV), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV2), dd.V.row(outV),
                            dd.H(inV2), dd.H(outV), dd.offset));
  dd.new_verts.push_back(split(dd.V.row(inV3), dd.V.row(outV),
                            dd.H(inV3), dd.H(outV), dd.offset));

  int h1 = newV_pos, h2 = newV_pos + 1, h3 = newV_pos + 2;
  // Create the new faces.
  std::vector<std::vector<int> > newFs;
  // Don't care about normals here... Can fix later (?)
  // Two simple faces (bottom and top).
  if (isValid({inV1, inV2, inV3}, dd.face_counts, dd.H)) {
    newFs.push_back({inV1, inV2, inV3});
    oTimes[0]++;
  }
  // TODO: here
  /*
  if ((inV1 == 270 || inV2 == 270 || inV3 == 270) &&
      (h1 == 2572 || h2 == 2572 || h3 == 2572))
    addOrig(T, 3, dd.faces, dd.faces_markers);
    */
  // Always add this one.
  newFs.push_back({h1, h2, h3});
  oTimes[1]++;
  // 6 additional faces.
  // Front: i1, h1, i2 + i2, h1, h2
  if (isValid({inV1, inV2, outV}, dd.face_counts, dd.H)) {
    newFs.push_back({inV1, newV_pos, inV2});
    newFs.push_back({inV2, newV_pos, newV_pos + 1});
    oTimes[2]++;
  }
  // Right: i2, i3, h2 + h2, i3, h3
  if (isValid({inV2, inV3, outV}, dd.face_counts, dd.H)) {
    newFs.push_back({inV2,         inV3, newV_pos + 1});
    newFs.push_back({newV_pos + 1, inV3, newV_pos + 2});
    oTimes[3]++;
  }
  // Left: i3, h3, i1 + i1, h3, h1
  if (isValid({inV3, inV1, outV}, dd.face_counts, dd.H)) {
    newFs.push_back({inV3, newV_pos + 2, inV1});
    newFs.push_back({inV1, newV_pos + 2, newV_pos});
    oTimes[4]++;
  }

  for (int i = 0; i < newFs.size(); ++i) {
    Eigen::RowVector3i t({newFs[i][0], newFs[i][1], newFs[i][2]});
    dd.faces.push_back(t);
    dd.faces_markers.push_back(3);
  }
}
}

// A tet is organized in the following fashion:
//
//       0
//       
//       3       (3 is in the background; 1 and 2 are in the foreground)
//  2        1
//
// So the faces (with counterclockwise normals) are:
// (0 1 3)
// (0 2 1)
// (3 2 0)
// (1 2 3)
//
//
// This method will perform three steps:
// 1. Add all faces, duplicating ones on the interior
// 2. Remove all duplicate verts (might have been caused during #1)
// 3. Remove all duplicate faces
void marching_tets(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  const Eigen::VectorXd& H,
  double offset,
  Eigen::MatrixXd& NV,
  Eigen::MatrixXi& NF,
  Eigen::VectorXi& I)
{
  using namespace Eigen;
  using namespace std;

  // Count the faces.
  std::map<std::vector<int>, int> face_counts;
  for (int i = 0; i < T.rows(); ++i) {
    std::vector<std::vector<int> > fs;
    fs.push_back({T(i, 0), T(i, 1), T(i, 3)});
    fs.push_back({T(i, 0), T(i, 2), T(i, 1)});
    fs.push_back({T(i, 3), T(i, 2), T(i, 0)});
    fs.push_back({T(i, 1), T(i, 2), T(i, 3)});
    for (auto &f : fs) {
      std::sort(f.begin(), f.end());
      // Add it to the map.
      face_counts[f]++;
    }
  }

  vector<Eigen::RowVector3i> faces;
  vector<int> faces_markers;
  vector<Eigen::RowVector3d> new_verts;
  int times[6];
  for (int i = 0; i < 6; i++) {
    times[i] = 0;
  }
  for (int i= 0; i < 5; ++i) {
    oTimes[i] = 0;
  }

  // Create data structure.
  MarchingTetsDat dd(V, H, faces, faces_markers, face_counts, new_verts, offset);

  // Check each tet face, add as needed.
  for (int i = 0; i < T.rows(); ++i) {
    // See if the tet is entirely inside.
    vector<int> inside, outside, identical, inside_t, outside_t;
    for (int j = 0; j < T.cols(); ++j) {
      if (H(T(i, j)) >= offset) {
        outside.push_back(j);
      } else if (H(T(i, j)) < offset) {
        inside.push_back(j);
      } else {
        identical.push_back(j);
      }

      if (H(T(i, j)) == GLOBAL::outside_temp) {
        outside_t.push_back(j);
      } else if (H(T(i, j)) == GLOBAL::inside_temp) {
        inside_t.push_back(j);
      }
    }

    // Ignore this tet if it's entirely outside.
    if (outside.size() == 4) continue;


    if (outside.size() == 0 && inside.size() == 0) {
      // degenerate, ignore.
      printf("WARNING: degenerate tet face found!!\n");
    } else if (outside.size() == 0) {
      times[0]++;
      // (A) Takes care of:
      //   inside: 1, identical: 3 (remove three duplicated faces later)
      //   inside: 2, identical: 2 (remove all four duplicated faces later)
      //   inside: 3, identical: 1 (remove all four duplicated faces later)
      //   inside: 4, identical: 0 (remove all four duplicated faces later)
      case0Out(dd, T.row(i), inside, outside, identical, inside_t, outside_t);
    } else if (identical.size() == 3 && outside.size() == 1) {
      times[1]++;
      // Don't do anything here. Will have already added it, so we don't want
      // to duplicate the face.
    } else if (inside.size() == 1) {
      // (these are colored green)
      times[2]++;
      // (B) Takes care of:
      //   inside: 1 outside: 3
      //   inside: 1 outside: 2 identical: 1
      //   inside: 1 outside: 1 identical: 2
      //
      case1In(dd, T.row(i), inside, outside, identical, inside_t, outside_t);
    } else if (inside.size() == 3 && outside.size() == 1) {
      // (these are colored orange)
      times[3]++;
      // (C) takes care of:
      //   inside: 3 outside: 1
      //
      case3In1Out(dd, T.row(i), inside, outside, identical, inside_t, outside_t);
    } else if (inside.size() == 2 && outside.size() > 0) {
      // (these are colored red)
      times[4]++;
      // (D) takes care of:
      //   inside: 2 outside: 1 identical: 1
      //   inside: 2 outside: 2 identical: 0
      //
      case2In2Out(dd, T.row(i), inside, outside, identical, inside_t, outside_t);
    } else {
      times[5]++;
      fprintf(stderr, "WARN: marching tets found something weird, with in:%lu out:%lu =:%lu\n",
              inside.size(), outside.size(), identical.size());
    }
  }

  printf("Finished marching tets with usages:\n");
  for (int i = 0; i < 6; ++i) {
    printf("  %d: %d\n", i, times[i]);
  }
  printf("Orange usages:\n");
  for (int i= 0; i < 5; ++i) {
    printf("  %d: %d\n", i, oTimes[i]);
  }

  // Copy verts
  NV.resize(V.rows() + new_verts.size(), 3);
  for (int i = 0; i < V.rows(); ++i) {
    NV.row(i) = V.row(i);
  }
  for (int i = 0; i < new_verts.size(); ++i) {
    NV.row(i + V.rows()) = new_verts[i];
  }
  // Set I
  I.resize(NV.rows());
  for (int i = 0; i < I.rows(); ++i) {
    if (i < V.rows()) {
      I(i) = i;
    } else {
      I(i) = -1;
    }
  }
  Eigen::VectorXi facesMarkers;
  facesMarkers.resize(faces.size());
  // Copy faces
  NF.resize(faces.size(), 3);
  for (int i = 0; i < faces.size(); ++i) {
    NF.row(i) = faces[i];
    facesMarkers(i) = faces_markers[i];
  }
  
  Eigen::MatrixXd newV;
  Eigen::MatrixXi newF;
  Eigen::VectorXi SVJ, SVI, I2;
  Helpers::viewTriMesh(NV, NF, facesMarkers);

  //igl::collapse_small_triangles(NV, NF, 1e-6, newF);
  //NF = newF;

  // Step 2: Remove all duplicate verts
  /*
  igl::remove_duplicate_vertices(NV, GLOBAL::EPS, newV, SVI, SVJ);
  I2.resize(newV.rows());
  for (int i = 0; i < newV.rows(); ++i) {
    I2(i) = I(SVI(i));
  }
  I = I2;
  for (int i = 0; i < NF.rows(); ++i) {
    for (int j = 0; j < NF.cols(); ++j) {
      NF(i, j) = SVJ(NF(i, j));
    }
  }
  NV = newV;
  //Helpers::viewTriMesh(NV, NF, I);
  */

  Helpers::removeDuplicates(NV, NF, I);

  igl::remove_unreferenced(NV, NF, newV, newF, SVI, SVJ);
  I2.resize(newV.rows());
  I2.setConstant(-1);
  for (int i = 0; i < I2.rows(); ++i) {
    I2(i) = I(SVJ(i));
  }
  I = I2;
  NV = newV;
  NF = newF;

  //Helpers::extractManifoldPatch(NV, NF, I);

  Helpers::viewTriMesh(NV, NF, I);


  // orient everything correctly.
  Eigen::VectorXi C;
  igl::orientable_patches(NF, C);
  igl::bfs_orient(NF, newF, C);
  NF = newF;
  igl::orient_outward(NV, NF, C, newF, SVJ);
  NF = newF;
}

