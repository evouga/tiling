#ifndef REMESH_H
#define REMESH_H

#include <Eigen/Core>

namespace Remeshing {
// Remesh the input mesh, attempting to get the target edge length e_length. If
// e_length is non-positive (including 0), then use the average edge length as
// the target.
void remesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &M,
            Eigen::MatrixXd &Vout, Eigen::MatrixXi &Fout, Eigen::VectorXi &Mout,
            double e_length = -1);
} // namespace Remeshing
#endif // REMESH_H
