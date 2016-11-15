#include "curvatureFlow.h"

#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/harmonic.h>
#include <igl/massmatrix.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"

static constexpr double gStoppingCriteria = 0.01;

#define FACE_DIM  3

void indices01(const Eigen::MatrixXi &F, const Eigen::VectorXi &orig,
               Eigen::VectorXi &I01) {

  // Temporary vector that will hold 1 if this is a 01 vertex (either on the
  // boundary or just one layer outside).
  Eigen::VectorXi temp_v(orig.rows());
  temp_v.setZero();
  // Find all the faces that are incident on orig.
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < FACE_DIM; ++j) {
      // If my marker is an original, everything next to me might be a +1
      // boundary vertex. Be over-zealous here, then adjust later.
      if (orig(F(i, j)) == GLOBAL::original_marker) {
        // Set all the other neighbors to 1
        for (int k = 1; k <= FACE_DIM - 1; ++k) {
          int idx = F(i, (j + k) % FACE_DIM);
          if (orig(idx) != GLOBAL::original_marker) {
            temp_v(idx) = 1;
          }
        }
      }
    }
  }

  // Also check for boundary markers (neighbors of +1 boundaries that are
  // also original vertices).
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < FACE_DIM; ++j) {
      // If I'm a +1 boundary vertex, then find all my neighbors that are also
      // original vertices.
      if (temp_v(F(i, j)) == 1) {
        for (int k = 1; k <= FACE_DIM - 1; k++) {
          int idx = F(i, (j + k) % FACE_DIM);
          // If I'm an original vertex, then I must be a boundary vertex.
          if (orig(idx) == GLOBAL::original_marker) {
            temp_v(idx) = 1;
          }
        }
      }
    }
  }
  
  // The number of 01 indices is the sum of the markers.
  int num_verts = temp_v.sum();
  I01.resize(num_verts);
  int temp_ctr = 0;
  // Create I01, which holds the indices of 01 vertices.
  for (int i = 0; i < temp_v.rows(); ++i) {
    if (temp_v(i)) {
      I01(temp_ctr++) = i;
    }
  }
}

void boundaryIndices(const Eigen::MatrixXi &F, const Eigen::VectorXi &orig,
                     Eigen::VectorXi &boundary) {
  Eigen::VectorXi temp_v(orig.rows());
  temp_v.setZero();

  for (int i = 0; i < F.rows(); ++i) {
    int num_inside = 0, num_outside = 0;

    // tri faces with 3 indices.
    for (int j = 0; j < FACE_DIM; ++j) {
      int idx = F(i, j);
      if (orig(idx) == GLOBAL::original_marker) num_inside++;
      else num_outside++;
    }
    
    // This face is both on the inside and the outside.
    if (num_inside && num_outside) {
      for (int j = 0; j < FACE_DIM; ++j) {
        int idx = F(i, j);
        if (orig(idx) == GLOBAL::original_marker) {
          temp_v(idx) = 1;
        }
      }
    }
  }

  // Create the final vector with just indices
  int num_v = temp_v.sum();
  boundary.resize(num_v);
  int temp_ctr = 0;
  for (int i = 0; i < temp_v.size(); ++i) {
    if (temp_v(i) == 1) {
      boundary(temp_ctr++) = i;
    }
  }
}

void outerIndices(const Eigen::VectorXi &orig,
                  Eigen::VectorXi &outside) {
  int num_v = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) != GLOBAL::original_marker) num_v++;
  }

  outside.resize(num_v);
  int outside_ctr = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) != GLOBAL::original_marker) {
      outside(outside_ctr++) = i;
    }
  }
}

void biharmonic_self(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                double timestep,
                Eigen::MatrixXd &Vg) {
  fprintf(stderr, "Made it inside this function!\n");
  // Get the indices of all the specific vertices.
  Eigen::VectorXi Itheta,  I0, I01, I, Ii;
  I = Eigen::VectorXi::LinSpaced(orig.rows(), 0, orig.rows() - 1);
  outerIndices(orig, Itheta);
  boundaryIndices(F, orig, I0);
  indices01(F, orig, I01);
  // theta_bar is just theta plus boundaries (I0)
  Eigen::VectorXi Itheta_bar(Itheta.rows() + I0.rows());
  Itheta_bar << Itheta, I0;
  // Ii is the complement of Itheta_bar
  Eigen::VectorXi IA; // Don't care about this, but needed for setdiff
  igl::setdiff(I,Itheta_bar, Ii, IA);

  fprintf(stderr,"Sizes are: Itheta:%d Itheta_bar:%d I0:%d I01:%d\n",
          Itheta.rows(), Itheta_bar.rows(), I0.rows(), I01.rows());
  // Get the entire Laplacian matrix.
  Eigen::SparseMatrix<double> L, M, Mdi;
  igl::cotmatrix(V, F, L);
  // Get the Laplacian matrices we care about.
  Eigen::SparseMatrix<double> Lttb, Ltbt, Ltb01;
  igl::slice(L, Itheta, Itheta_bar, Lttb);
  igl::slice(L, Itheta_bar, Itheta, Ltbt);
  igl::slice(L, Itheta_bar, I01, Ltb01);

  // Updated vertices (all of them)
  Eigen::MatrixXd Vc = V;
  // Vertices that don't change
  Eigen::MatrixXd ui;
  igl::slice(V, Ii, 1, ui);
  // Only known vertices
  Eigen::MatrixXd uf01;
  igl::slice(V, I01, 1, uf01);
  // Only unknown vertices.
  Eigen::MatrixXd u;


  // Set up index_into, which will tell the vector and index for each vertex.
  Eigen::MatrixXi index_into(V.rows(), 2);
  for (int i = 0; i < Ii.rows(); ++i) {
    index_into(Ii(i), 0) = 0;
    index_into(Ii(i), 1) = i;
  }
  for (int i = 0; i < I01.rows(); ++i) {
    index_into(I01(i), 0) = 1;
    index_into(I01(i), 1) = i;
  }
  for (int i = 0; i < Itheta.rows(); ++i) {
    index_into(Itheta(i), 0) = 2;
    index_into(Itheta(i), 1) = i;
  }
  
  igl::viewer::Viewer v;
  v.callback_key_down = [&](igl::viewer::Viewer& v, unsigned char key, int modifier) {
    // Press spacebar to get the next version.
    if (key == ' ') {

      // The voronoi hybrid lumped mass matrix is better than the barycentric
      // version.
      igl::massmatrix(Vc, F, igl::MASSMATRIX_TYPE_VORONOI, M);
      // Only care about entries of unknown vertices
      igl::slice(M, Itheta_bar, Itheta_bar, Mdi);
      // Want the inverse of this.
      for (int i = 0; i < Mdi.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Mdi, i); it; ++it) {
          it.valueRef() = 1.0 / it.value();
        }
      }

      const auto &LHS = Lttb * Mdi * Ltbt;
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(LHS);
      u = solver.solve(-Lttb * Mdi * Ltb01 * uf01);

      assert(solver.info() == Eigen::Success);
      // Calculate the difference.
      double diff = 0;
      // Vertices updated like:
      for (int i = 0; i < index_into.rows(); ++i) {
        int vec = index_into(i, 0);
        int idx = index_into(i, 1);
        if (vec == 0) {
          diff += (Vc.row(i) - ui.row(idx)).norm();
          Vc.row(i) = ui.row(idx);
        } else if (vec == 1) {
          diff += (Vc.row(i) - uf01.row(idx)).norm();
          Vc.row(i) = uf01.row(idx);
        } else {
          diff += (Vc.row(i) - u.row(idx)).norm();
          Vc.row(i) = u.row(idx);
        }
      }

      printf("Difference this round is %f\n", diff);

      v.data.set_vertices(Vc);
    }
    return true;
  };

  v.data.set_mesh(Vc, F);
  v.launch();

}

void convertVectorToIndices(const Eigen::VectorXi &orig,
                            Eigen::VectorXi &b) {
  // b will contain all the boundary vertices.
  int num_orig = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) == GLOBAL::original_marker) {
      num_orig++;
    }
  }
  b.resize(num_orig);
  int orig_i = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) == GLOBAL::original_marker) {
      b(orig_i++) = i;
    }
  }
}

void biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                Eigen::MatrixXd &Vc) {
  Eigen::VectorXi b;
  // Get the correct indices.
  convertVectorToIndices(orig, b);
  // The positions of these indices are held constant (just keep the input values)
  Eigen::MatrixXd V_bc = igl::slice(V, b, 1);

  // Initialize Vc with input vertices V
  Vc = V;

  igl::viewer::Viewer v;
  v.callback_key_down = [&](igl::viewer::Viewer& v, unsigned char key, int modifier) {
    // Press spacebar to get the next version.
    if (key == ' ') {
      // How much should we change by?
      Eigen::MatrixXd D;
      Eigen::MatrixXd V_bc = igl::slice(Vc, b, 1);
      igl::harmonic(Vc,F, b,V_bc, 2, D);

      // Calculate the difference.
      double diff = 0;
      // Vertices updated like:
      for (int i = 0; i < V.rows(); ++i) {
        diff += (Vc.row(i) - D.row(i)).norm();
      }
      printf("Difference this round is %f\n", diff);
      
      //Vc = Vc + D;
      Vc = D;

      v.data.set_vertices(Vc);
    }
    return true;
  };

  printf("\nViewing biharmonic mesh.\n");
  printf("  Press ' ' (space) to deform shape progressively\n");
  printf("  Close down image to return shape to previous function\n");
  v.data.set_mesh(Vc, F);
  v.launch();

}


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

  // Set the original rows to zero in the Laplace matrix
  L.prune([&orig](int r, int c, float) {
    return (orig(r) != GLOBAL::original_marker);// || (r == c);
  });

  igl::viewer::Viewer v;
  do {
    // Compute the mass matrix
    igl::massmatrix(Vc, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
    //igl::massmatrix(Vc, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    //for (int i = 0; i < V.rows(); ++i) {
    //  M.coeffRef(i,i) = 0;
    //}
    printf("Number of nonzeros for M: %d vs rows %d\n", M.nonZeros(), M.rows());
    // Solve (M-delta*L) U = M*U
    const auto &S = (M - timestep*L);
    //Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver(S);
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver(S);
    U = solver.solve(M*Vc).eval();
    assert(solver.info() == Eigen::Success);

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
      difference += (U.row(i) - Vc.row(i)).norm();
    }
    printf("difference this round is %lf (stop is %lf)\n",
           difference, gStoppingCriteria);
    Vc = U;

    v.data.set_mesh(Vc, F);
    v.launch();

  } while (difference > gStoppingCriteria);
}
