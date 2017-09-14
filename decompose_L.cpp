#include "decompose_L.h"

#include <vector>

#include <igl/cotmatrix_entries.h>

/**
 * Decompose L into D^T\star D.
 *
 * @param V - vertices of mesh
 * @param F - faces of mesh
 */
void decompose_L(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 Eigen::SparseMatrix<double> &D, Eigen::SparseMatrix<double> &star) {
  // Construct D
  D.resize(F.rows() * 3, V.rows());
  //                                 # non-zeros per vertex
  D.reserve(Eigen::VectorXi::Constant(2, V.rows()));
  //D.reserve(F.rows()*3 * 2);
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int from = F(i, j);
      int to = F(i, (j + 1) % 3);
      D.insert(i*3 + j, from) = -1;
      D.insert(i*3 + j, to) = 1;
    }
  }

  // Construct star from cot_entries.
  Eigen::MatrixXd cot_entries;
  // cot_entries is #Fx3 where each row corresponds to edge:
  //   ( [1,2], [2,0], [0,1] )
  igl::cotmatrix_entries(V,F, cot_entries);

  star.resize(F.rows() * 3, F.rows() * 3);
  star.reserve(Eigen::VectorXi::Constant(F.rows() * 3, 1));
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int e = i * 3 + j;
      // Need the negative value for this to work.
      star.coeffRef(e, e) -= cot_entries(i, (j + 2) % 3);
    }
  }

  //fprintf(stderr, "D is %ld,%ld and star is %ld,%ld\nFinished with matrix\n",
  //        D.rows(), D.cols(), star.rows(), star.cols());
}
