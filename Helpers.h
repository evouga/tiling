#ifndef HELPERS_H
#define HELPERS_H

#include <Eigen/Core>

namespace Helpers {

void removeDuplicates(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O);

Eigen::VectorXi getFaceMarkers(const Eigen::MatrixXi& F, const Eigen::VectorXi& M);

// TODO(bradyz): Create a tri-mesh data structure.
void tetrahedralize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    const Eigen::VectorXi &VM, const Eigen::VectorXi &FM,
                    Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                    Eigen::MatrixXi &TF, Eigen::VectorXi &TO);

void triangulate(const Eigen::MatrixXd &P, const Eigen::MatrixXi &E,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F);

} // namespace Helpers

#endif // HELPERS_H
