#include "curvatureFlow.h"

#include <iostream> // cout

#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/harmonic.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"

static constexpr double gStoppingCriteria = 0.0;

#define FACE_DIM  3

namespace {
// Extracts a vector with only the indices of original vertices.
void convertVectorToIndices(const Eigen::VectorXi &orig,
                            Eigen::VectorXi &b) {
  // b will contain all the boundary vertices.
  int num_orig = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) != GLOBAL::nonoriginal_marker) {
      num_orig++;
    }
  }
  b.resize(num_orig);
  int orig_i = 0;
  for (int i = 0; i < orig.rows(); ++i) {
    if (orig(i) != GLOBAL::nonoriginal_marker) {
      b(orig_i++) = i;
    }
  }
}

// Changes all the faces in the interior (i.e. those that have no non-original
// neighbors) to non-original faces. Maintains a single layer of fixed
// vertices.
void removeInterior(const Eigen::MatrixXi &F, const Eigen::VectorXi &orig,
                    Eigen::VectorXi &newOrig) {
  typedef Eigen::Array<bool,Eigen::Dynamic,1> VectorXb;
  VectorXb nonorig_neighbor(orig.rows());
  nonorig_neighbor.setConstant(false);

  newOrig = orig;

  // Look at all the rows
  for (int i = 0; i < F.rows(); ++i) {
    bool contains_no = false;
    // For each face, check and see if there are any nonoriginal markers
    for (int j = 0; j < F.cols(); ++j) {
      int vi = F(i, j);
      if (orig(vi) == GLOBAL::nonoriginal_marker) {
        contains_no = true;
        break;
      }
    }

    // if the face does contain a non-original, mark all the vertices correspondingly.
    if (contains_no) {
      for (int j = 0; j < F.cols(); ++j) {
        int vi = F(i, j);
        nonorig_neighbor(vi) = true;
      }
    } else {
      for (int j = 0; j < F.cols(); ++j) {
        int vi = F(i, j);
      }
    }
  }

  for (int i = 0; i < newOrig.rows(); ++i) {
    if (!nonorig_neighbor(i)) {
      newOrig(i) = GLOBAL::nonoriginal_marker;
    }
  }
}
} // namespace

double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &orig,
                  Eigen::MatrixXd &Vc, double change_val,
                  bool remove_interior) {
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);

  return biharmonic(V, F, orig, L, Vc, change_val, remove_interior);
}

double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &orig, const Eigen::SparseMatrix<double> &L,
                  Eigen::MatrixXd &Vc, double change_val,
                  bool remove_interior) {
  // First thing we need to do: remove interior vertices
  Eigen::VectorXi new_orig = orig;
  if (remove_interior) {
    removeInterior(F, orig, new_orig);
  }
  Eigen::VectorXi b;
  // Get the correct indices.
  convertVectorToIndices(new_orig, b);
  // The positions of these indices are held constant (just keep the input values)
  Eigen::MatrixXd V_bc = igl::slice(V, b, 1);

  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  // Initialize Vc with input vertices V
  Vc = V;

  for (int counter = 0; counter < 10; counter++) {
    // Recompute M, leave L alone.
    igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

    // How much should we change by?
    Eigen::MatrixXd D;
    igl::harmonic(L,M, b,V_bc, 2, D);

    // Calculate the difference in vertex positions.
    double diff = 0;
    // Vertices updated like:
    for (int i = 0; i < V.rows(); ++i) {
      diff += (Vc.row(i) - D.row(i)).norm();
    }

    // Update the values.
    Vc = D;
  }

  // Calculate the energy.
  Eigen::SparseMatrix<double> Mi;
  igl::invert_diag(M, Mi);
  return (Vc.transpose() * L * Mi * L * Vc).trace();
}

// Previous function that allows one to view the biharmonic.
void biharmonic_view(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                Eigen::MatrixXd &Vc, bool remove_interior) {
  Eigen::VectorXi new_orig = orig;
  if (remove_interior) {
    removeInterior(F, orig, new_orig);
  }

  Eigen::VectorXi b;
  // Get the correct indices.
  convertVectorToIndices(new_orig, b);
  // The positions of these indices are held constant (just keep the input values)
  Eigen::MatrixXd V_bc = igl::slice(V, b, 1);

  Eigen::SparseMatrix<double> L, M;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  // Initialize Vc with input vertices V
  Vc = V;

  igl::viewer::Viewer v;
  v.callback_key_down = [&](igl::viewer::Viewer& v, unsigned char key, int modifier) {
    if (key == 'L') {
      // Recompute the Laplacian
      igl::cotmatrix(Vc,F,L);
      printf("Recomputing Laplacian...\n");
    } else if (key == 'R') {
      printf("Resetting..\n");
      Vc = V;
      v.data.set_vertices(V);
      return true;
    } else if (key == ' ') {
      // Press spacebar to get the next version.
      
      // Recompute M, leave L alone.
      igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
      //igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_VORONOI,M);
      //igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);

      // How much should we change by?
      Eigen::MatrixXd D;
      //Eigen::MatrixXd V_bc = igl::slice(Vc, b, 1);
      igl::harmonic(L,M, b,V_bc, 2, D);
      //igl::harmonic(Vc,F, b,V_bc, 2, D);

      // Calculate the difference.
      double diff = 0;
      // Vertices updated like:
      for (int i = 0; i < V.rows(); ++i) {
        diff += (Vc.row(i) - D.row(i)).norm();
      }
      printf("Difference this round is %f\n", diff);
      
      // Update the values.
      //Vc = Vc + D;
      Vc = D;

      // Calculate the energy.
      Eigen::SparseMatrix<double> Mi;
      igl::invert_diag(M,Mi);
      auto en = (Vc.transpose() * L * Mi * L * Vc).trace();
      std::cout << "Energy is " << en << std::endl;
      
      // Set the colors.
      Eigen::MatrixXd cols;
      igl::jet(new_orig, true, cols);
      v.data.set_vertices(Vc);
      v.data.set_colors(cols);
    }
    return true;
  };

  printf("Mesh has %ld verts and %ld faces\n", Vc.rows(), F.rows());
  printf("\nViewing biharmonic mesh.\n");
  printf("  Press ' ' (space) to deform shape progressively\n");
  printf("  Press 'L' to (re)compute the Laplacian\n");
  printf("  Press 'R' to reset the vertices\n");
  printf("  Close down image to return shape to previous function\n");
  v.data.set_mesh(Vc, F);
  Eigen::MatrixXd cols;
  igl::jet(new_orig, true, cols);
  v.data.set_vertices(Vc);
  v.data.set_colors(cols);
  v.launch();
}


// include the vertices, faces, and markers. Output new vertices.
void computeCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &orig,
                          double timestep,
                          Eigen::MatrixXd &Vc) {
  // Compute the gradient/laplacian
  Eigen::SparseMatrix<double> L, Lo, M;
  igl::cotmatrix(V, F, L);
  Eigen::MatrixXd U;
  Vc = V;
  double difference;
  
  Lo = L;
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
      
    // Update the values.
    Vc = U;

    // Calculate the energy.
    Eigen::SparseMatrix<double> Mi;
    igl::invert_diag(M,Mi);
    auto en = (Vc.transpose() * Lo * Mi * Lo * Vc).trace();
    std::cout << "Energy is " << en << std::endl;

    v.data.set_mesh(Vc, F);
    v.launch();

  } while (difference > gStoppingCriteria);
}
