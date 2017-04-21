// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Nathan Clement <nathanlclement@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// Will perform marching tets on a tet-mesh, and return the offset surface
// at a given amount, delta.
// 
// Inputs:
//   V: #V by 3 list of vertices
//   T: #T by 4 list of tet indices into V
//   H: #V by 1 list of values at each vertex
//   offset: the offset surface level
//
// Outputs:
//   NV: #NV by 3 list of new vertices
//   NF: #NF by 3 list of face indices into NV
//   I: #NV by 3 list of indices such that if NV(i) != -1, NV(i) == V(I(i))
void marching_tets(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  const Eigen::VectorXd& H,
  double offset,
  Eigen::MatrixXd& NV,
  Eigen::MatrixXi& NF,
  Eigen::VectorXi& I);
