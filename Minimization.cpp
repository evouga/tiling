#include "Minimization.h"

#include <ctime>

#include <dlib/optimization.h>

#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_face_normals.h>

#include "curvatureFlow.h"

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
typedef dlib::matrix<double,0,1> column_vector;
typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;

namespace {

class test_function_deriv {
 public:
  test_function_deriv(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const std::vector<int> *to_ignore) 
      : _Vo(V), _Fo(F), _to_ignore(to_ignore)
  {
    // Set up the movable vector.
    _movable.resize(V.rows());
    _movable.setConstant(true);
    _num_movable = V.rows();

    if (to_ignore) {
      // if it's defined, set each position to false.
      for (int idx : *to_ignore) {
        _movable(idx) = false;
        // Also keep track of this.
        _num_movable--;
      }
    }
  }

  // This will return the derivative of the energy function, given the input
  // values.
  //const column_vector deriv(const column_vector& arg) const {
  const column_vector operator()(const column_vector& arg) const {
    // Compute these matrices first.
    Eigen::MatrixXd V = _Vo;
    // Need to get this special (full) matrix.
    int mov_idx = 0;
    for (int i = 0; i < V.rows(); ++i) {
      if (_movable(i)) {
        for (int j = 0; j < 3; ++j) {
          V(i, j) = arg(mov_idx * 3 + j);
        }
        mov_idx++;
      } else {
        for (int j = 0; j < 3; ++j) {
          V(i, j) = _Vo(i, j);
        }
      }
    }
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V,_Fo, L);
    igl::massmatrix(V,_Fo, igl::MASSMATRIX_TYPE_DEFAULT, M);
    for (int i = 0 ; i < V.rows(); ++i) {
      if (!_movable(i)) 
        M.coeffRef(i, i) = std::numeric_limits<double>::infinity();
    }
    igl::invert_diag(M, Minv);

    // Get each derivative separately, then add them all up.
    Eigen::MatrixXd dv = deriv_dL(V,L,Minv) + deriv_dV(V,L,Minv) 
        + deriv_dM(V,_Fo,L,Minv) + deriv_dTheta(V,L,Minv);


    // Convert them back into single-column format.
    mov_idx = 0;
    column_vector res(arg.size());
    for (int i = 0; i < V.rows(); ++i) {
      if (!_movable(i)) continue;
      for (int j = 0; j < 3; ++j) {
        res(mov_idx * 3 + j) = dv(i, j);
      }
      mov_idx++;
    }

    return res;
  }

 private:
  const Eigen::MatrixXd deriv_dL(
      const Eigen::MatrixXd &V,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    // TODO
    Eigen::MatrixXd ret(V.rows(), 3);
    ret.setZero();
    return ret;
  }

  // The contribution to the gradient $\nabla_{v_k} E$ is
  // $$4\left(L^{\intercal} M^{-1}L V\right)_k^T.$$
  const Eigen::MatrixXd deriv_dV(
      const Eigen::MatrixXd &V,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    return  4 * (L.transpose() * Mi * L * V).transpose();
  }

  // The contribution from the area to the gradient $\nabla_{v_k} E$ is
  // thus:
  //    $$\frac{1}{3}\sum_{f\supset v_k} 
  //          \left(\sum_{v_i\subset f} \|(M^{-1}LV)_i\|^2\right) 
  //                \left(v^f_+-v^f_-\right)\times\mathbf{n}_f.$$
  const Eigen::MatrixXd deriv_dM(
      const Eigen::MatrixXd &V,
      const Eigen::MatrixXi &F,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    // Need to create mapping of vertex->face
    std::vector<std::vector<int> > v2face;
    v2face.resize(V.rows());

    for (int i = 0; i < F.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        v2face[F(i, j)].push_back(i);
      }
    }

    // Pre-compute this matrix.
    auto milv = Mi * L * V;

    // Also precompute face normals.
    Eigen::MatrixXd n;
    igl::per_face_normals(V,F, n);

    // Finally, compute the value.
    Eigen::MatrixXd ret(V.rows(), 3);
    ret.setZero(); // null it out.
    for (int i = 0; i < V.rows(); ++i) {
      // Look at each face that contains i
      for (int fi : v2face[i]) {
        const auto &f = F.row(fi);
        // Compute the multiplier from each of the three vertices.
        double mult = 0;
        // Look at each vertex that is in f
        for (int j = 0; j < 3; ++j) {
          mult += milv.row(f(j)).squaredNorm();
        }

        // Then, look at the other two vertices attached to this face.
        // Compute next and prev vertices.
        int next, prev;
        if (i == f(0)) {
          prev = f(2); next = f(1); 
        } else if (i == f(1)) {
          prev = f(0); next = f(2); 
        } else { // i == f(2)
          prev = f(1); next = f(0); 
        }

        // Finally, compute the cross product with the face normal.
        const Eigen::RowVector3d diff = V.row(next) - V.row(prev);
        const Eigen::RowVector3d f_n = n.row(fi);
        Eigen::RowVector3d dv = diff.cross(f_n);
        ret.row(i) += mult * dv;
      }
    }

    return ret * -1. / 3;
  }

  const Eigen::MatrixXd deriv_dTheta(
      const Eigen::MatrixXd &V,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    // TODO
    Eigen::MatrixXd ret(V.rows(), 3);
    ret.setZero();
    return ret;
  }

  const Eigen::MatrixXd &_Vo; // Won't ever change.
  const Eigen::MatrixXi &_Fo; // Won't ever change.
  const std::vector<int> *_to_ignore; // Won't ever change.
  VectorXb _movable; // Won't ever change. Just indices.
  int _num_movable; // Number of mutable vertices.
};

class test_function {
 public:
  test_function(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const std::vector<int> *to_ignore) 
      : _Vo(V), _Fo(F), _to_ignore(to_ignore)
  {
    // Set up the movable vector.
    _movable.resize(V.rows());
    _movable.setConstant(true);
    _num_movable = V.rows();

    if (to_ignore) {
      // if it's defined, set each position to false.
      for (int idx : *to_ignore) {
        _movable(idx) = false;
        // Also keep track of this.
        _num_movable--;
      }
    }
  }

  column_vector getStartingPoint() const {
    column_vector start(_num_movable * 3);
    int new_idx = 0;;
    
    for (int i = 0; i < _Vo.rows(); ++i) {
      if (!_movable(i)) continue;

      // Update all three vertices.
      for (int j = 0; j < 3; ++j) {
        start(new_idx * 3 + j) = _Vo(i, j);
      }
      // Then use the next index.
      new_idx++;
    }
    return start;
  }

  Eigen::MatrixXd getFinalVerts(const column_vector& arg) const {
    int new_idx = 0;
    Eigen::MatrixXd Vt;
    Vt.resize(_Vo.rows(), _Vo.cols());
    for (int i = 0; i < _Vo.rows(); ++i) {
      // Ignore stationary things.
      if (!_movable(i)) {
        for (int j = 0; j < 3; ++j) {
          Vt(i, j) = _Vo(i, j);
        }
        continue;
      }

      for (int j = 0; j < 3; ++j) {
        Vt(i, j) = arg(new_idx * 3 + j);
      }
      new_idx++;
    }

    return Vt;
  }

  double operator() ( const column_vector& arg) const {
    std::clock_t here_start = std::clock();

    Eigen::MatrixXd Vt = getFinalVerts(arg);
    double en = biharmonic_energy(Vt, _Fo, _to_ignore);
    //printf("Inside function, value is %lf\n", en);

    // Return the biharmonic energy of the two.
    //double en = biharmonic_energy(Vt, _Fo, _to_ignore);
    return en;
  }

  void print(const column_vector& arg) const {
    Eigen::MatrixXd Vt = getFinalVerts(arg);
    int new_idx = 0;
    for (int i = 0; i < _Vo.rows(); ++i) {
      if (fabs((Vt.row(i) - _Vo.row(i)).norm()) < 1e-6) continue;

      printf("%d: ", i);
      for (int j = 0; j < 3; ++j) {
        if (_movable(i)) {
          printf("%lf", arg(new_idx * 3 + j));
        } else {
          printf("na");
        }
        printf(",%lf,%lf(%lf) ", Vt(i, j), _Vo(i, j),
               Vt(i, j) - _Vo(i, j));
      }
      // Then increment the new index if we need to.
      if (_movable(i)) new_idx++;
      printf("\n");
    }
  }

 private:
  const Eigen::MatrixXd &_Vo; // Won't ever change.
  const Eigen::MatrixXi &_Fo; // Won't ever change.
  const std::vector<int> *_to_ignore; // Won't ever change.
  VectorXb _movable; // Won't ever change. Just indices.
  int _num_movable; // Number of mutable vertices.
};

} // namespace

double minimize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const std::vector<int> *to_ignore /* can be null */,
                Eigen::MatrixXd &minV) {

  printf("Starting minimize function...\n");
  test_function tf(V, F, to_ignore);
  test_function_deriv tf_dv(V, F, to_ignore);
  printf("Just initialized tf function\n");
  column_vector start = tf.getStartingPoint();
  printf("Created starting point, energy is %lf\n", tf(start));
  tf.print(start);
  printf("Number of dimensions is %ld\n", start.size());
  auto approx_deriv = dlib::derivative(tf)(start);
  auto ours = tf_dv(start);
  std::cout << "Difference between analytic derivative and numerical approximation of derivative: "
            << length(approx_deriv - ours) << std::endl;
  std::cout << approx_deriv << "\nAnd ours:\n" << ours << std::endl;

  find_min_using_approximate_derivatives(
      //dlib::lbfgs_search_strategy(100),
      dlib::bfgs_search_strategy(),
      dlib::objective_delta_stop_strategy(1e-1).be_verbose(),
      tf, start, -1);
  printf("Finished approximate derivatives\n");

  minV = tf.getFinalVerts(start);
  printf("got final value (energy is %lf\n", tf(start));

  return tf(start);
}
