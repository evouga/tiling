#ifndef DECOMPOSE_L_H
#define DECOMPOSE_L_H

#include <Eigen/Core>
#include <Eigen/Sparse>

/*
 * Will decompose L into it's Hodge-star definition, which is:
 *
 *   L = D^T \star D
 *
 * where L is the |V|x|V| cotangent Laplacian, D is |V|x|E|, and
 * \star is |E|x|E|.
 */
void decompose_L(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 Eigen::SparseMatrix<double> &D, Eigen::SparseMatrix<double> &star);
#endif // DECOMPOSE_L_H
