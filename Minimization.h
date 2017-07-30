#ifndef MINIMIZATION_H
#define MINIMIZATION_H

#include <vector>

#include <Eigen/Core>

double minimize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const std::vector<int> *to_ignore /* can be null */,
                Eigen::MatrixXd &minV);

#endif // MINIMIZATION_H
