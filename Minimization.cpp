#include "Minimization.h"

#include <ctime>

#include <dlib/optimization.h>

#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_face_normals.h>

#include "curvatureFlow.h"
#include "decompose_L.h"

// In dlib, the general purpose solvers optimize functions that take a column
// vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
typedef dlib::matrix<double,0,1> column_vector;
typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
typedef Eigen::Matrix<double, -1, -1, Eigen::RowMajor> RowMatrixXd;
typedef Eigen::Matrix<int, -1, -1, Eigen::RowMajor> RowMatrixXi;

static int GLOB_ROW_i;

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
    
    // Need to create mapping of vertex->face
    _v2face.resize(V.rows());

    for (int i = 0; i < F.rows(); ++i) {
      for (int j = 0; j < 3; ++j) {
        _v2face[F(i, j)].push_back(i);
      }
    }

    // Construct _uE and _EMAP
    Eigen::MatrixXi E, uE2E; // Don't care about these.
    //igl::unique_edge_map(F, E, _uE, _EMAP, uE2E);
  }

  // This will return the derivative of the energy function, given the input
  // values.
  //const column_vector deriv(const column_vector& arg) const {
  const column_vector operator()(const column_vector& arg) const {
    // Compute these matrices first; make sure it's in row-major order!!
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> V = _Vo;
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
    igl::massmatrix(V,_Fo, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    for (int i = 0 ; i < V.rows(); ++i) {
      if (!_movable(i)) 
        M.coeffRef(i, i) = std::numeric_limits<double>::infinity();
    }
    igl::invert_diag(M, Minv);

    // Get each derivative separately, then add them all up.
    auto dL = deriv_dL(V,_Fo,L,Minv),
         dV = deriv_dV(V,L,Minv),
         dM = deriv_dM(V,_Fo,L,Minv),
         dTh= deriv_dTheta(V,L,Minv);
    Eigen::MatrixXd dv = dL + dV + dM + dTh;
    /*
    std::cout << "dL: " << dL.colwise().minCoeff()
                   << " max " << dL.colwise().maxCoeff() << std::endl;
    std::cout << "dV: " << dV.colwise().minCoeff()
                   << " max " << dV.colwise().maxCoeff() << std::endl;
    std::cout << "dM: " << dM.colwise().minCoeff()
                   << " max " << dM.colwise().maxCoeff() << std::endl;
    std::cout << "dTh: " << dTh.colwise().minCoeff()
                   << " max " << dTh.colwise().maxCoeff() << std::endl;
                   */
    // TESTING TESTING just set to dM
    //Eigen::MatrixXd dv = solo_dM(V,_Fo,L,Minv);

    /*
    fprintf(stderr, "Just about to decompose L...\n");
    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> star;
    decompose_L(V,_Fo, D,star);
    // VxE * ExE * ExV
    Eigen::MatrixXd decL = D.transpose() * star * D;

    int N = dv.rows();
    Eigen::MatrixXd together(N, 5*3);
    together.block(0,0, N,3) = deriv_dL(V,L,Minv);
    together.block(0,3, N,3) = deriv_dV(V,L,Minv);
    together.block(0,6, N,3) = deriv_dM(V,_Fo,L,Minv);
    together.block(0,9, N,3) = deriv_dTheta(V,L,Minv);
    together.block(0,12, N,3) = dv;
    std::cout << "All derivs are:\n"
         << together << std::endl;
         */

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
  // Will get the next index and previous vertex index, and the zero-based
  // index of the two in f.
  void getNextPrev(const Eigen::RowVector3i &f, int idx,
                   int &next, int &prev, int &sp, int &sm, int &s0,
                   bool swap = false) const {
    if (idx == f(0)) {
      prev = f(2); next = f(1); 
      sm = 2; sp = 1; s0 = 0;
    } else if (idx == f(1)) {
      prev = f(0); next = f(2); 
      sm = 0; sp = 2; s0 = 1;
    } else { // idx == f(2)
      prev = f(1); next = f(0); 
      sm = 1; sp = 0; s0 = 2;
    }
    // Need to do this to match derivatives of dM (why??)
    if (swap) {
      int temp = prev;
      prev = next; next = temp;
      temp = sm;
      sm = sp;
      sp = temp;
    }
  }
  // Wrapper function.
  void getNextPrev(const Eigen::RowVector3i &f, int idx,
                   int &next, int &prev) const {
    int sp, sm, s0;
    getNextPrev(f,idx, next,prev, sp,sm,s0, false);

  }

  const Eigen::MatrixXd deriv_dL(
      const Eigen::MatrixXd &V,
      const Eigen::MatrixXi &F,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    Eigen::MatrixXd ret(V.rows(), 3);
    ret.setZero();

    // Also precompute face normals and area.
    Eigen::MatrixXd n, F_2area;
    igl::per_face_normals(V,F, n);
    igl::doublearea(V,F, F_2area);
    n = -n;

    // Also decompose L into D^T\starD
    Eigen::SparseMatrix<double> D, star;
    decompose_L(V,F, D,star);
    
    // And pre-compute some other matrices.
    const RowMatrixXd DMiLV = D * Mi * L * V;
    const RowMatrixXd DV = D * V;

    // Then compute function
    for (int i = 0; i < V.rows(); ++i) {
      //fprintf(stderr, "Working on row %d\n", i);
      // Three components for each face.
      for (int fi : _v2face[i]) {
        //fprintf(stderr, "Working on f %d\n", fi);
        const auto &f = F.row(fi);
        const Eigen::Vector3d &fn = n.row(fi);
        double Af = 0.5 * F_2area(fi);  // A(f) = 0.5 * 2*A(f)
        double Af2 = Af * Af * 2.;      // 2 * A(f)^2
        //fprintf(stderr, "Area is %lf\n", Af);

        // Get next and prev indices.
        int next, prev, sm, sp, s0;
        getNextPrev(f, i, next, prev, sm, sp, s0);
        int spkf, smkf, s0kf;
        if (i == f(0)) {
          spkf = 0; smkf = 2; s0kf = 1; next = f(1); prev = f(2);
        } else if (i == f(1)) {
          spkf = 1; smkf = 0; s0kf = 2; next = f(2); prev = f(0);
        } else { // if (i == f(2))
          spkf = 2; smkf = 1; s0kf = 0; next = f(0); prev = f(1);
        }
        /*
        // Need to reverse this.
        if (i == f(0)) {
          spkf = 2; smkf = 0; s0kf = 1; next = f(2); prev = f(1);
        } else if (i == f(1)) {
          spkf = 0; smkf = 1; s0kf = 2; next = f(0); prev = f(2);
        } else { // if (i == f(2))
          spkf = 1; smkf = 2; s0kf = 0; next = f(1); prev = f(0);
        }
        */

        // First component:
        double c1 = DMiLV.row(fi * 3 + spkf).dot(
                              DV.row(fi * 3 + spkf)); // s_plus
        // (v_+ - v_-) / A(f)
        auto v1_1 = (V.row(next) - V.row(prev)) / Af;
        //// (v_- - v_k) / A(f)
        //auto v1_1 = (V.row(prev) - V.row(i)) / Af;
        // (v_k - v_-) * (v_+ - v_-)
        // ------------------------- (v_+ - v_-)
        //         2 A(f)^2
        Eigen::Vector3d v1_2a = 
            ((V.row(i) - V.row(prev)) * 
             (V.row(next) - V.row(prev)) / Af2) *
            (V.row(next) - V.row(prev));
        // Need to do some funky shit to get this to have correct dimensions.
        Eigen::RowVector3d v1_2 = (v1_2a.transpose().cross(fn)).transpose();
        Eigen::RowVector3d v1 = c1 * (v1_1 + v1_2);

        // Second term
        double c2 = DMiLV.row(fi * 3 + smkf).dot(
                              DV.row(fi * 3 + smkf)); // s_minus
        // (v_- - v_+) / A(f)
        auto v2_1 = (V.row(prev) - V.row(next)) / Af;
        // (v_+ - v_k) / A(f)
        //auto v2_1 = (V.row(next) - V.row(i)) / Af;
        // (v_k - v_+) * (v_- - v_+)
        // ------------------------- (v_+ - v_-)
        //         2 A(f)^2
        Eigen::Vector3d v2_2a = 
            ((V.row(i) - V.row(next)) * 
             (V.row(prev) - V.row(next)) / Af2) *
            (V.row(next) - V.row(prev));
        Eigen::RowVector3d v2_2 = (v2_2a.transpose().cross(fn)).transpose();
        Eigen::RowVector3d v2 = c2 * (v2_1 + v2_2);

        // Third term.
        double c3 = DMiLV.row(fi * 3 + s0kf).dot(
                              DV.row(fi * 3 + s0kf)); // s_original
        // (2v_k - v_+ - v_-) / A(f)
        auto v3_1 = (2 * V.row(i) - V.row(next) - V.row(prev)) / Af;
        // (v_+ - v_k) * (v_- - v_k)
        // ------------------------- (v_+ - v_-)
        //         2 A(f)^2
        Eigen::Vector3d v3_2a = 
            ((V.row(next) - V.row(i)) * 
             (V.row(prev) - V.row(i)) / Af2) *
            (V.row(next) - V.row(prev));
        Eigen::RowVector3d v3_2 = (v2_2a.transpose().cross(fn)).transpose();
        Eigen::RowVector3d v3 = c3 * (v3_1 + v3_2);


        // Total is sum of parts -- summed over all faces.
        ret.row(i) += v1 + v2 + v3;
      }
    }

    return ret;
  }

  // The contribution to the gradient $\nabla_{v_k} E$ is
  // $$4\left(L^{\intercal} M^{-1}L V\right)_k^T.$$
  const Eigen::MatrixXd deriv_dV(
      const Eigen::MatrixXd &V,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    return  4.0 * (L.transpose() * Mi * L * V);
  }

  Eigen::MatrixXd mat_cp(const Eigen::RowVector3d &v) const {
    Eigen::MatrixXd v_x(3, 3);
    v_x << 0, -v(2), v(1),
           v(2), 0, -v(0),
           -v(1), v(0), 0;
    return v_x;
  }

  // Derivative of mass matrix (sanity check).
  const Eigen::MatrixXd solo_dM(
      const Eigen::MatrixXd &V,
      const Eigen::MatrixXi &F,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi) const {
    Eigen::MatrixXd ret(V.rows(), 3);
    ret.setZero();
    
    // Also precompute face normals.
    Eigen::MatrixXd n;
    igl::per_face_normals(V,F, n);
    n = -n;

    // Then compute function.
    for (int i = 0; i < V.rows(); ++i) {
      Eigen::Vector3d dMii;
      dMii << 0, 0, 0;
      // Look at each face that contains i
      for (int fi : _v2face[i]) {
        const auto &f = F.row(fi);

        // Then, look at the other two vertices attached to this face.
        // Compute next and prev vertices.
        int next, prev;
        getNextPrev(f, i, next, prev);

        // Then compute the multiplier.
        //dMii += n.row(fi).transpose() * 
        dMii += n.row(fi) *
            ( mat_cp(V.row(i) - V.row(prev)) - 
              mat_cp(V.row(i) - V.row(next)) );
      }

      ret.row(i) = dMii.transpose();
    }

    return 1. / 6. * ret;

  }

  // The contribution from the area to the gradient $\nabla_{v_k} E$ is
  // thus:
  //    $$\frac{1}{3}\sum_{f\supset v_k} 
  //          \left(\sum_{v_i\subset f} \|(M^{-1}LV)_i\|^2\right) 
  //                \left(v^f_+-v^f_-\right)\times\mathbf{n}_f.$$
  const Eigen::MatrixXd deriv_dM(
      const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &V,
      const Eigen::MatrixXi &F,
      const Eigen::SparseMatrix<double> &L,
      const Eigen::SparseMatrix<double> &Mi, bool solo_dv = false) const {
    // Pre-compute this matrix.
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> milv = Mi * L * V;

    // Also precompute face normals.
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> n;
    igl::per_face_normals(V,F, n);

    // Finally, compute the value.
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> ret(V.rows(), 3);
    ret.setZero(); // null it out.
    for (int i = 0; i < V.rows(); ++i) {
      // Look at each face that contains i
      for (int fi : _v2face[i]) {
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
        getNextPrev(f, i, next,prev);

        // Finally, compute the cross product with the face normal.
        const Eigen::Vector3d diff = (V.row(next) - V.row(prev));
        const Eigen::Vector3d f_n = n.row(fi);
        //Eigen::RowVector3d dv = diff.cross(f_n);
        ret.row(i) += mult * diff.cross(f_n);
        //ret(i, 0) += mult * 
            //((V.row(next) - V.row(prev)).transpose().cross(
                    //n.row(fi).transpose()));
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
  std::vector<std::vector<int> > _v2face; // Won't ever change.
  Eigen::MatrixXi _uE; // List of unique undirected edges.
  Eigen::VectorXi _EMAP; // #F*3 mapping from F to _uE
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

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(Vt,_Fo, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    // TESTING!!
    ////double ret =  M.diagonal().sum();
    //return M.coeffRef(GLOB_ROW_i, GLOB_ROW_i);


    Eigen::SparseMatrix<double> Mi;
    igl::invert_diag(M, Mi);
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(Vt, _Fo, L);
    double en = (Vt.transpose() * L.transpose() * Mi * L * Vt).trace();

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
  /*
  auto ours = tf_dv(start);
  GLOB_ROW_i = 0; // This one is bad, before.
  std::cout << " working on " << GLOB_ROW_i << std::endl;
  auto approx_deriv = dlib::derivative(tf)(start);
  for (int i = GLOB_ROW_i+1; i < V.rows(); ++i) {
    bool same = true;
    for (int j = 0; j < 3; ++j) {
      if (std::abs(approx_deriv((i-1) * 3 + j) -
                    ours((i-1) * 3 + j)) > 1e-6) {
        same = false;
      }
      printf("%lf ", approx_deriv((i-1) * 3 + j));
    }
    printf("vs");
    for (int j = 0; j < 3; ++j) {
      printf(" %lf", ours((i-1) * 3 + j));
    }
    printf("%s\n", same ? " :)" : "***");
    std::cout << " working on " << i << std::endl;
    GLOB_ROW_i = i;
    start = tf.getStartingPoint();
    approx_deriv = dlib::derivative(tf)(start);
  }
  std::cout << "Difference between analytic derivative and numerical approximation of derivative: "
            << length(approx_deriv - ours) << std::endl;
  for (int i = 0; i < approx_deriv.size(); ++i) {
    printf("%c%9.6lf %9.6lf\n", i % 3 == 0 ? '*' : ' ', approx_deriv(i), ours(i));
  }
  //std::cout << approx_deriv << "\nAnd ours:\n" << ours << std::endl;
  */
  /*
  find_min_using_approximate_derivatives(
      //dlib::lbfgs_search_strategy(100),
      dlib::bfgs_search_strategy(),
      dlib::objective_delta_stop_strategy(1e-3).be_verbose(),
      tf, start, -1);
      */
  /* This uses our gradient. */
  /*
  find_min(dlib::bfgs_search_strategy(),
           dlib::objective_delta_stop_strategy(1e-1).be_verbose(),
           tf, tf_dv, start, -1);
  auto approx_deriv = dlib::derivative(tf)(start);
  auto ours = tf_dv(start);
  for (int i = 0; i < approx_deriv.size(); ++i) {
    printf("%c%9.6lf %9.6lf\n", i % 3 == 0 ? '*' : ' ', approx_deriv(i), ours(i));
  }
           */
  double en = tf(start);
  double en2 = en;
  auto dv = tf_dv(start);
  auto prev = start;
  double start_en = en;
  double delta = 1e-8;
  int num = 0;
  while (length(dv) > 0 && num < 200) {// && en <= en2) {
    prev = start;
    // Move in a small amount that direction.
    start -= dv * delta;
    // Recompute the energy.
    en2 = en;
    en = tf(start);
    // Compute the derivatives
    dv = tf_dv(start);
    printf("[%d] Energy is %lf (change %lf; diff is %lf)\n", 
           num++, en, en - start_en, length(dv));
  }
  if (en2 < en) {
    printf("Switching back to prev\n");
    start = prev;
  }

  /*
  auto approx_deriv_final = dlib::derivative(tf)(start);
  auto ours_final = tf_dv(start);
  for (int i = 0; i < approx_deriv_final.size(); ++i) {
    printf("%c%9.6lf %9.6lf\n", i % 3 == 0 ? '*' : ' ', 
           approx_deriv_final(i), ours_final(i));
  }
  */

  printf("Finished approximate derivatives\n");

  minV = tf.getFinalVerts(start);
  printf("got final value (energy is %lf, vert diffs to follow)\n", tf(start));
  return tf(start);
}
