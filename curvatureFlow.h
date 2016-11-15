#ifndef CURVATURE_FLOW_H
#define CURVATURE_FLOW_H

#include <Eigen/Core>

void computeCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &M,
                          double timestep,
                          Eigen::MatrixXd &Vc);
// Given some input mesh V,F, produce an output Vc
void biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &M,
                          Eigen::MatrixXd &Vc);
#endif // CURVATURE_FLOW_H
