#ifndef CURVATURE_FLOW_H
#define CURVATURE_FLOW_H

#include <Eigen/Core>

void computeCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &M,
                          double timestep,
                          Eigen::MatrixXd &Vc);

// Given some input mesh V,F, produce an output Vc
// Also includes a flag to remove any interior "original" markers, so the
// normals are not messed up.
void biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                Eigen::MatrixXd &Vc, double change_val = 10,
                bool remove_interior = true);
// Same function as above, for viewing.
void biharmonic_view(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                Eigen::MatrixXd &Vc,
                bool remove_interior = true);

#endif // CURVATURE_FLOW_H
