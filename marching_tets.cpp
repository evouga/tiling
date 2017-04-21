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
  for (int i = 0; i < 3; ++i) {
    newP(i) = (how_much * a(i) + (1 - how_much) * b(i));
  }

  return newP;
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
  std::vector<int> t(3);
  for (int i = 0; i < 3; ++i) {
    t[i] = f(i);
  }
  return getMapCounts(counts, t);
}
bool isValid(std::vector<int> &c,
             const std::map<std::vector<int>, int> &counts,
             const Eigen::VectorXd &H) {
  /*
  for (int i : c) {
    if (H(i) == GLOBAL::outside_temp) return false;
  }
  */

  return getMapCounts(counts, c) == 1;
}
bool isValid(const Eigen::RowVector3i &f,
             const std::map<std::vector<int>, int> &counts,
             const Eigen::VectorXd &H) {
  for (int i = 0; i < 3; ++i) {
    //if (H(f(i)) == GLOBAL::outside_temp) return false;
  }

  return getMapCounts(counts, f) == 1;
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
  map<std::vector<int>, int> face_counts;
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
  int oTimes[5];
  for (int i= 0; i < 5; ++i) {
    oTimes[i] = 0;
  }

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

    // A couple of cases to be concerned about with temperature.
    if (inside_t.size() > 0 && outside_t.size() > 0 && 0) {
      switch(outside_t.size()) {
        case 4:
          // degenerate.
          printf("WARNING: degenerate tet face found (all same temperature)!!\n");
          continue;
        case 1: {
          times[1]++;
          // If only one is on the boundary, make a face out of the other three.
          Eigen::RowVector3i newF;
          int fc = 0;
          for (int j = 0; j < 4; ++j) {
            if (j == outside_t[0]) continue;
            newF(fc++) = T(i, j);
          }
          faces.push_back(newF);
          faces_markers.push_back(1); // not used.
          // Don't do anything else.
          continue;
        }
        case 2:
        case 3:
          // Don't do anything to this face (collapse both verts).
          continue;
      } 
    }

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
      // Add all 4 faces.
      Eigen::RowVector3i f1, f2, f3, f4;
      f1 << T(i, 0), T(i, 1), T(i, 3);
      f2 << T(i, 0), T(i, 2), T(i, 1);
      f3 << T(i, 3), T(i, 2), T(i, 0);
      f4 << T(i, 1), T(i, 2), T(i, 3);
      if (isValid(f1, face_counts, H)) {
        faces.push_back(f1);
        faces_markers.push_back(0);
      }
      if (isValid(f2, face_counts, H)) {
        faces.push_back(f2);
        faces_markers.push_back(0);
      }
      if (isValid(f3, face_counts, H)) {
        faces.push_back(f3);
        faces_markers.push_back(0);
      }
      if (isValid(f4, face_counts, H)) {
        faces.push_back(f4);
        faces_markers.push_back(0);
      }
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
      // Only one inside, need to divide each of the remaining three in half,
      // then add one new face (ignore the other three)
      int inV = T(i, inside[0]);
      int outV1, outV2, outV3;
      if (outside.size() == 3) {
        outV1 = T(i, outside[0]);
        outV2 = T(i, outside[1]);
        outV3 = T(i, outside[2]);
      } else if (outside.size() == 2) {
        outV1 = T(i, outside[0]);
        outV2 = T(i, outside[1]);
        outV3 = T(i, identical[0]);
      } else if (outside.size() == 1) {
        outV1 = T(i, outside[0]);
        outV2 = T(i, identical[0]);
        outV3 = T(i, identical[1]);
      }

      int newV_pos = V.rows() + new_verts.size();
      new_verts.push_back(split(V.row(inV), V.row(outV1),
                                H(inV), H(outV1), offset));
      new_verts.push_back(split(V.row(inV), V.row(outV2),
                                H(inV), H(outV2), offset));
      new_verts.push_back(split(V.row(inV), V.row(outV3),
                                H(inV), H(outV3), offset));
      
      // Create the new faces.
      Eigen::RowVector3i newF1, newF2, newF3, newF4;
      // Don't care about normals here... Can fix later (?)
      newF1 << newV_pos, newV_pos + 1, newV_pos + 2;
      newF2 << inV, newV_pos, newV_pos + 1;
      newF3 << inV, newV_pos + 1, newV_pos + 2;
      newF4 << inV, newV_pos + 2, newV_pos;
      // This face is automatically added.
      faces.push_back(newF1);
      faces_markers.push_back(2);
      // Only add the rest of the faces if they're on a a boundary.
      if (isValid({inV, outV1, outV2}, face_counts, H)) {
        faces.push_back(newF2);
        faces_markers.push_back(2);
      }
      if (isValid({inV, outV2, outV3}, face_counts, H)) {
        faces.push_back(newF3);
        faces_markers.push_back(2);
      }
      if (isValid({inV, outV3, outV1}, face_counts, H)) {
        faces.push_back(newF4);
        faces_markers.push_back(2);
      }
    } else if (inside.size() == 3 && outside.size() == 1) {
      // (these are colored orange)
      times[3]++;
      // (C) takes care of:
      //   inside: 3 outside: 1
      //
      // Only one outside, need to divide each of the remaining three in half,
      // then add one new face.
      int outV = T(i, outside[0]);
      int inV1 = T(i, inside[0]);
      int inV2 = T(i, inside[1]);
      int inV3 = T(i, inside[2]);

      int newV_pos = V.rows() + new_verts.size();
      new_verts.push_back(split(V.row(inV1), V.row(outV),
                                H(inV1), H(outV), offset));
      new_verts.push_back(split(V.row(inV2), V.row(outV),
                                H(inV2), H(outV), offset));
      new_verts.push_back(split(V.row(inV3), V.row(outV),
                                H(inV3), H(outV), offset));
      
      // Create the new faces.
      std::vector<std::vector<int> > newFs;
      // Don't care about normals here... Can fix later (?)
      // Two simple faces (bottom and top).
      if (isValid({inV1, inV2, inV3}, face_counts, H)) {
        newFs.push_back({inV1, inV2, inV3});
        oTimes[0]++;
      }
      // Always add this one.
      newFs.push_back({newV_pos, newV_pos + 1, newV_pos + 2});
      oTimes[1]++;
      // 6 additional faces.
      // Front: i1, h1, i2 + i2, h1, h2
      if (isValid({inV1, inV2, outV}, face_counts, H)) {
        newFs.push_back({inV1, newV_pos, inV2});
        newFs.push_back({inV2, newV_pos, newV_pos + 1});
        oTimes[2]++;
      }
      // Right: i2, i3, h2 + h2, i3, h3
      if (isValid({inV2, inV3, outV}, face_counts, H)) {
        newFs.push_back({inV2,         inV3, newV_pos + 1});
        newFs.push_back({newV_pos + 1, inV3, newV_pos + 2});
        oTimes[3]++;
      }
      // Left: i3, h3, i1 + i1, h3, h1
      if (isValid({inV3, inV1, outV}, face_counts, H)) {
        newFs.push_back({inV3, newV_pos + 2, inV1});
        newFs.push_back({inV1, newV_pos + 2, newV_pos});
        oTimes[4]++;
      }

      for (int i = 0; i < newFs.size(); ++i) {
        Eigen::RowVector3i t;
        t << newFs[i][0], newFs[i][1], newFs[i][2];
        faces.push_back(t);
        faces_markers.push_back(3);
      }
    } else if (inside.size() == 2 && outside.size() > 0) {
      // (these are colored red)
      times[4]++;
      // (D) takes care of:
      //   inside: 2 outside: 1 identical: 1
      //   inside: 2 outside: 2 identical: 0
      //
      // 2 inside, 2 outside. Insert four faces between 4 points.
      int outV1 = T(i, outside[0]);
      int outV2 = T(i, outside[1]);
      int inV1 = T(i, inside[0]);
      int inV2 = T(i, inside[1]);

      // Create 4 new points.
      int newV_pos = V.rows() + new_verts.size();
      new_verts.push_back(split(V.row(inV1), V.row(outV1), // h11
                                H(inV1), H(outV1), offset));
      new_verts.push_back(split(V.row(inV2), V.row(outV1), // h21
                                H(inV2), H(outV1), offset));
      new_verts.push_back(split(V.row(inV1), V.row(outV2), // h12
                                H(inV1), H(outV2), offset));
      new_verts.push_back(split(V.row(inV2), V.row(outV2), // h22
                                H(inV2), H(outV2), offset));

      // Create two faces
      std::vector<std::vector<int> > newFs;

      // Front: i2, h22, h21
      if (isValid({inV2, outV2, outV1}, face_counts, H)) {
        newFs.push_back({inV2, newV_pos + 3, newV_pos + 1});
      }
      // Back: i1, h11, h12
      if (isValid({inV1, outV1, outV2}, face_counts, H)) {
        newFs.push_back({inV1, newV_pos, newV_pos + 2});
      }
      // Side: i2, i1, h22 + h12, i1, h22
      if (isValid({inV1, inV2, outV2}, face_counts, H)) {
        newFs.push_back({inV2,         inV1, newV_pos + 3});
        newFs.push_back({newV_pos + 3, inV1, newV_pos + 2});
      }
      // Side: i2, i1, h21 + h21, i1, h11
      if (isValid({inV1, inV2, outV1}, face_counts, H)) {
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
        Eigen::RowVector3i t;
        t << newFs[i][0], newFs[i][1], newFs[i][2];
        faces.push_back(t);
        faces_markers.push_back(4);
      }
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
  Helpers::viewTriMesh(NV, NF, I);
  */
  Helpers::viewTriMesh(NV, NF, facesMarkers);

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
  Helpers::viewTriMesh(NV, NF, I);

  //Helpers::extractManifoldPatch(NV, NF, I);

  // orient everything correctly.
  Eigen::VectorXi C;
  igl::orientable_patches(NF, C);
  igl::bfs_orient(NF, newF, C);
  NF = newF;
  igl::orient_outward(NV, NF, C, newF, SVJ);
  NF = newF;
}
