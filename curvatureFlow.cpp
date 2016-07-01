#include "curvatureFlow.h"

#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"

static constexpr double gStoppingCriteria = 0.01;

// include the vertices, faces, and markers. Output new vertices.
void computeCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &orig,
                          double timestep,
                          Eigen::MatrixXd &Vc) {
  // Compute the gradient/laplacian
  Eigen::SparseMatrix<double> L, M;
  igl::cotmatrix(V, F, L);
  Eigen::MatrixXd U;
  Vc = V;
  double difference;

  do {
    // Compute the mass matrix
    igl::massmatrix(Vc, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    // Solve (M-delta*L) U = M*U
    const auto &S = (M - timestep*L);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(S);
    assert(solver.info() == Eigen::Success);
    U = solver.solve(M*Vc).eval();

    /*
    // Compute centroid and subtract
    Eigen::VectorXd dblA;
    igl::doublearea(U,F,dblA);

    double area = 0.5*dblA.sum();
    Eigen::MatrixXd BC;
    igl::barycenter(U,F,BC);
    Eigen::RowVector3d centroid(0,0,0);
    for(int i = 0;i<BC.rows();i++) {
      centroid += 0.5*dblA(i)/area*BC.row(i);
    }
    U.rowwise() -= centroid;

    // Normalize to unit surface area (important for numerics)
    U.array() /= sqrt(area);
    */
    // Calculate difference between U and Vc
    difference = 0;
    for (int i = 0; i < U.rows(); ++i) {
      if (orig(i) == GLOBAL::original_marker) {
        // Keep this the same.
        U.row(i) = V.row(i);
      } else {
        difference += (U.row(i) - Vc.row(i)).norm();
      }
    }
    printf("difference this round is %lf (stop is %lf)\n",
           difference, gStoppingCriteria);
    Vc = U;
    igl::viewer::Viewer v;
    v.data.set_mesh(Vc, F);
    v.launch();

  } while (difference > gStoppingCriteria);
}
