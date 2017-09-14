#ifndef DECOMPOSE_L_H
#define DECOMPOSE_L_H

#include <Eigen/Core>
#include <Eigen/Sparse>

/*
 * Will decompose L into its Hodge-star definition, which is:
 *
 *   L = D^T \star D
 *
 * where L is the |V|x|V| cotangent Laplacian, D is |V|x|E|, and
 * \star is |E|x|E|. Note that each edge in E is directed (i.e.
 * |E| = |F| * 3, whereas for undirected edges, |E| != |F|*3.
 
 * @param V - vertices of mesh
 * @param F - faces of mesh
 * @param OUT D, star
 */
void decompose_L(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 Eigen::SparseMatrix<double> &D, Eigen::SparseMatrix<double> &star);

#endif // DECOMPOSE_L_H
