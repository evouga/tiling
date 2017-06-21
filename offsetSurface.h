#ifndef OFFSET_SURFACE_H
#define OFFSET_SURFACE_H

#include <Eigen/Core>

#include <igl/AABB.h>

namespace OffsetSurface {

class Triangulation {
 public:
  igl::AABB<Eigen::MatrixXd,3> _tree; // AABB tree
  Eigen::MatrixXd _V; // Vertices
  Eigen::MatrixXi _T; // Tets
  Eigen::VectorXd _C; // Vertex values
};

void generateOffsetSurface(
    Triangulation &T,
    double off, Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff);

void generateOffsetSurface(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &TT, const Eigen::VectorXd &C,
    double off, Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff,
    Triangulation &T
    );

// Returns offset vertices, faces, and original markers.
/*
void generateOffsetSurface_naive(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO,
    const Eigen::VectorXd &C, double off,
    Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff, Eigen::VectorXi &Ooff);
    */

// Offset surface using marching tets.
void marchingOffsetSurface(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO,
    const Eigen::VectorXd &C, double off,
    Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff, Eigen::VectorXi &Ooff);
} // namespace OffsetSurface

#endif // OFFSET_SURFACE_H
