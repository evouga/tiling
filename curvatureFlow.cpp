#include "curvatureFlow.h"

#include <iostream>
#include <set>
#include <map>

#include <igl/adjacency_list.h>
#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/harmonic.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/remove_unreferenced.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "Helpers.h"

static constexpr double gStoppingCriteria = 0.0;

#define FACE_DIM  3

using namespace std;

namespace {

// Extracts a vector with only the indices of original vertices.
Eigen::VectorXi extract_fixed_vertices(const Eigen::VectorXi &O) {
  // b will contain all the boundary vertices.
  int number_fixed = 0;
  for (int i = 0; i < O.rows(); ++i) {
    if (O(i) != GLOBAL::nonoriginal_marker)
      number_fixed++;
  }

  Eigen::VectorXi fixed_vertices = Eigen::VectorXi(number_fixed);

  int offset = 0;
  for (int i = 0; i < O.rows(); ++i) {
    if (O(i) != GLOBAL::nonoriginal_marker)
      fixed_vertices(offset++) = i;
  }

  return fixed_vertices;
}

// Changes all the faces in the interior (i.e. those that have no non-original
// neighbors) to non-original faces.
// Maintains a single layer of fixed vertices.
Eigen::VectorXi remove_interior(const Eigen::MatrixXi &F,
                                const Eigen::VectorXi &O) {
  // Initialize.
  Eigen::VectorXi O_new = O;

  vector<bool> not_original(O.rows(), false);
  for (int i = 0; i < F.rows(); ++i) {
    bool has_nonoriginal_neighbor = false;

    // For each face, check and see if there are any nonoriginal markers.
    for (int j = 0; j < F.cols(); ++j) {
      if (O(F(i, j)) == GLOBAL::nonoriginal_marker)
        has_nonoriginal_neighbor = true;
    }

    // Change the marker.
    if (has_nonoriginal_neighbor) {
      for (int j = 0; j < F.cols(); ++j)
        not_original[F(i, j)] = true;
    }
  }

  for (int i = 0; i < O_new.rows(); ++i) {
    if (!not_original[i])
      O_new(i) = GLOBAL::nonoriginal_marker;
  }

  return O_new;
}

void triangulate_boundary(const Eigen::MatrixXd &V,
                          const vector<vector<int> > &graph,
                          map<int, map<int, int> > &edge_count,
                          set<int> &visited,
                          int start_node, int offset, bool extend_top,
                          vector<Eigen::VectorXi> &F_border) {
  // Make sure to set the start as visited.
  visited.insert(start_node);

  // A continous path of border vertices' indices.
  vector<int> path;

  vector<int> to_visit;
  to_visit.push_back(start_node);

  while (to_visit.size() > 0) {
    int u = to_visit.back();
            to_visit.pop_back();
    path.push_back(u);

    for (int v : graph[u]) {
      if (visited.find(v) != visited.end())
        continue;
      else if (edge_count[u][v] != 1)
        continue;

      visited.insert(v);
      
      // Make sure to break to ensure continuous path.
      to_visit.push_back(v);
      break;
    }
  }

  // The new faces that connect the old and new vertices.
  vector<Eigen::VectorXi> F_side;
  for (int i = 0; i < path.size(); ++i) {
    int j = (i + 1) % path.size();

    // Create two new faces. Normals are not correct.
    Eigen::Vector3i f1(offset + path[i], offset + path[j], path[i]);
    Eigen::Vector3i f2(path[j], path[i], offset + path[j]);

    F_side.push_back(f1);
    F_side.push_back(f2);
  }

  // Figure out whether the path went clockwise or counterclockwise;
  int clockwise = 0;
  int counter_clockwise = 0;
  for (int i = 0; i < path.size(); ++i) {
    const Eigen::Vector3d &a = V.row(path[i]);
    const Eigen::Vector3d &b = V.row(path[(i+1) % path.size()]);
    const Eigen::Vector3d &c = V.row(path[(i+2) % path.size()]);

    Eigen::Vector3d u = (c - b);
    Eigen::Vector3d v = (a - b);
    Eigen::Vector3d w = u.cross(v);

    if (w(2) > GLOBAL::EPS)
      clockwise++;
    else if (w(2) < GLOBAL::EPS)
      clockwise--;
  }

  // Flip face normals if not oriented correctly.
  if ((clockwise > counter_clockwise) == extend_top) {
    for (Eigen::VectorXi &face : F_side) {
      int tmp = face(0);
      face(0) = face(1);
      face(1) = tmp;
    }
  }

  // Add the faces for this connected component.
  for (Eigen::VectorXi &face : F_side)
    F_border.push_back(face);
}

void extendVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    const Eigen::VectorXi &O,
                    Eigen::MatrixXd &V_new, Eigen::MatrixXi &F_new,
                    Eigen::VectorXi &O_new,
                    bool extend_top,
                    double ext_amt=0.1) {
  double boundary_z;
  if (extend_top)
    boundary_z = V.colwise().maxCoeff()(2);
  else
    boundary_z = V.colwise().minCoeff()(2);

  // Identify extreme vertices
  set<int> V_boundary_outer;
  set<int> V_boundary_inner;
  set<int> top_v_all;

  // Go through vertices and look for points on boundary.
  for (int i = 0; i < V.rows(); ++i) {
    double z = V(i, 2);

    if (fabs(z - boundary_z) < GLOBAL::EPS) {
      if (O(i) != GLOBAL::nonoriginal_marker)
        V_boundary_outer.insert(i);
      else
        V_boundary_inner.insert(i);

      top_v_all.insert(i);
    }
  }

  // Degenerate case when there is only one maximal point.
  if (V_boundary_outer.size() < 1) {
    V_new = V;
    F_new = F;
    O_new = O;
    return;
  }

  // edge_count[u][v] = 1 means this is an edge on the boundary.
  map<int, map<int, int> > edge_count;

  // Use boundary vertices to find all the faces that lie on the boundary.
  for (int i = 0; i < F.rows(); i++) {
    bool all_vertices_on_boundary = true;
    for (int j = 0; j < 3; j++) {
      if (top_v_all.find(F(i, j)) == top_v_all.end())
        all_vertices_on_boundary = false;
    }

    if (!all_vertices_on_boundary)
      continue;

    // Construct undirected graph.
    for (int j = 0; j < 3; j++) {
      int u = F(i, j);
      int v = F(i, (j+1) % 3);

      edge_count[u][v]++;
      edge_count[v][u]++;
    }
  }

  // Used to offset into raised points.
  int V_size = V.rows();

  // Duplicate and extend all points (can use speedup).
  Eigen::MatrixXd V_with_duplicates(V.rows() * 2, V.cols());

  // Vertex i is duplicated and extended at V.rows() + i.
  for (int i = 0; i < V.rows(); ++i) {
    V_with_duplicates.row(i) = V.row(i);
    V_with_duplicates.row(V_size + i) = V.row(i);

    if (extend_top)
      V_with_duplicates(V_size + i, 2) += ext_amt;
    else
      V_with_duplicates(V_size + i, 2) -= ext_amt;
  }

  // Need to triangulate the sides by traversing the border.
  vector<Eigen::VectorXi> F_border;
  set<int> visited;
  vector<vector<int> > graph;
  igl::adjacency_list(F, graph);

  // DFS along the boundary for all connected components.
  for (int node : V_boundary_outer) {
    if (visited.find(node) != visited.end())
      continue;

    triangulate_boundary(V, graph, edge_count, visited, node, V_size,
                         extend_top, F_border);
  }

  // Populate F_new. Includes faces of raised points and new triangulated sides.
  F_new.resize(F.rows() + F_border.size(), F.cols());

  // Look at the original points to see if the vertices should be raised.
  for (int i = 0; i < F.rows(); ++i) {
    // Copy the original face.
    F_new.row(i) = F.row(i);

    bool should_raise = false;

    // All inner faces are shifted.
    for (int j = 0; j < F.cols(); ++j) {
      if (V_boundary_inner.find(F(i, j)) != V_boundary_inner.end())
        should_raise = true;
    }

    // Also allowed if all 3 vertices are on the boundary.
    if (!should_raise &&
        top_v_all.find(F(i, 0)) != top_v_all.end() &&
        top_v_all.find(F(i, 1)) != top_v_all.end() &&
        top_v_all.find(F(i, 2)) != top_v_all.end()) {
      should_raise = true;
    }

    // Use the vertices that have been offset.
    if (should_raise) {
      for (int j = 0; j < 3; ++j)
        F_new(i, j) += V_size;
    }
  }

  // Add the additional faces.
  for (int i = 0; i < F_border.size(); ++i)
    F_new.row(F.rows() + i) = F_border[i];

  // J is new indices into F, I new indices into V
  Eigen::VectorXi I, J; 
  igl::remove_unreferenced(V_with_duplicates.rows(), F_new, J, I);

  // Update faces after removing vertices.
  for (int i = 0; i < F_new.rows(); ++i) {
    for (int j = 0; j < 3; ++j)
      F_new(i, j) = J(F_new(i, j));
  }

  // Update V, with only things from I.
  igl::slice(V_with_duplicates, I, 1, V_new);

  // Update the markers.
  O_new.resize(V_new.rows());
  for (int i = 0; i < O_new.rows(); ++i) {
    int updated_index = I(i);

    // Get the original marker.
    if (updated_index >= V.rows())
      O_new(i) = O(updated_index - V_size);
    else
      O_new(i) = O(updated_index);
  }
}

// Wrapper for IGL function to provide fixed values for each dimension
// separately.
bool igl_harmonic_single(const Eigen::SparseMatrix<double> &Q,
                         const Eigen::VectorXi &b, const Eigen::VectorXd &V,
                         int n,
                         Eigen::VectorXd &W) {
  W.resize(n);

  // Minimize trace( 0.5*W' * Q * W + constant)
  // subject to:
  //    W(b:) = V, and
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(),
                                      true, data);

  const Eigen::VectorXd B = Eigen::VectorXd::Zero(n);

  return min_quad_with_fixed_solve(data, B, V, Eigen::VectorXd(), W);
}

double biharmonic_energy(const Eigen::MatrixXd &V,
                         const Eigen::MatrixXi &F) {
  // Laplacian.
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Mass matrix.
  Eigen::SparseMatrix<double> M;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

  Eigen::SparseMatrix<double> Mi;
  igl::invert_diag(M, Mi);

  return (V.transpose() * L.transpose() * Mi * L * V).trace();
}

} // namespace

double biharmonic_new(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const Eigen::VectorXi &O,
                      Eigen::MatrixXd &V_new, Eigen::MatrixXi &F_new,
                      Eigen::VectorXi &O_new) {
  // The updated markers only keep the boundary vertices.
  Eigen::VectorXi O_exterior = remove_interior(F, O);

  Eigen::MatrixXd V_prepared;
  Eigen::MatrixXd V_tmp;
  Eigen::MatrixXi F_tmp;
  Eigen::VectorXi O_tmp;

  // Remember which vertices are the max/min so we can track the added ones.
  double z_max_prev = V.colwise().maxCoeff()(2);
  double z_min_prev = V.colwise().minCoeff()(2);

  // Add vertices to top and bottom (respectively).
  extendVertices(V, F, O_exterior,
                 V_tmp, F_tmp, O_tmp, true);
  extendVertices(V_tmp, F_tmp, O_tmp,
                 V_prepared, F_new, O_new, false);

  // Mark the ones that are out of the original scale (new vertices)
  Eigen::VectorXi xy_nonfixed = O_new;

  for (int i = 0; i < V_prepared.rows(); ++i) {
    double z = V_prepared(i, 2);

    // If this is a new vertex, allow the x, y to move.
    if (z > z_max_prev || z < z_min_prev)
      xy_nonfixed(i) = GLOBAL::nonoriginal_marker;
  }

  // Also get the Laplacian
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V_prepared, F_new, L);

  // No need to remove interior.
  return biharmonic(V_prepared, F_new, O_new, xy_nonfixed, L,
                    V_new);
}

double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &O,
                  Eigen::MatrixXd &V_new, bool remove_interior_vertices) {
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  Eigen::VectorXi O_new = O;

  if (remove_interior_vertices)
    O_new = remove_interior(F, O);

  return biharmonic(V, F, O_new, O_new, L, V_new);
}

double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &O, const Eigen::VectorXi &fixed_xy,
                  const Eigen::SparseMatrix<double> &L,
                  Eigen::MatrixXd &V_new) {
  // Get the indices that are fixed.
  Eigen::VectorXi b_xy = extract_fixed_vertices(fixed_xy);
  Eigen::VectorXi b_z = extract_fixed_vertices(O);

  // Initialize V_new with input vertices V.
  V_new = V;

  // Mass matrix, harmonic weights.
  Eigen::SparseMatrix<double> M;
  Eigen::SparseMatrix<double> Q;

  for (int counter = 0; counter < 10; counter++) {
    // The positions of these indices are held constant (keep the input values).
    Eigen::MatrixXd V_bc = igl::slice(V_new, b_xy, 1);
    Eigen::VectorXd V_x_fixed = V_bc.col(0);
    Eigen::VectorXd V_y_fixed = V_bc.col(1);
    Eigen::VectorXd V_z_fixed = igl::slice(V_new, b_z, 1).col(2);

    // Recompute M based on the new vertex positions, leave L alone.
    igl::massmatrix(V_new, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    // How much should we change by?
    /* Previous version, all change the same.
    Eigen::MatrixXd D;
    igl::harmonic(L,M, b,V_bc, 2, D);
    Eigen::VectorXd Dx = D.col(0), Dy = D.col(1), Dz = D.col(2);
    */

    // Get biharmonic (k=2) weights.
    igl::harmonic(L, M, 2, Q);

    // Solve for each of these.
    Eigen::VectorXd V_x_new, V_y_new, V_z_new;
    igl_harmonic_single(Q, b_xy, V_x_fixed, V.rows(), V_x_new);
    igl_harmonic_single(Q, b_xy, V_y_fixed, V.rows(), V_y_new);
    igl_harmonic_single(Q,  b_z, V_z_fixed, V.rows(), V_z_new);

    // Calculate the difference in vertex positions.
    double diff = 0.0;
    for (int i = 0; i < V.rows(); ++i) {
      Eigen::RowVector3d new_position(V_x_new(i), V_y_new(i), V_z_new(i));
      diff += (V_new.row(i) - new_position).norm();
    }

    // Update the values.
    V_new.col(0) = V_x_new;
    V_new.col(1) = V_y_new;
    V_new.col(2) = V_z_new;
  }

  // Calculate the energy.
  return biharmonic_energy(V_new, F);
}

// Previous function that allows one to view the biharmonic.
void biharmonic_view(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXi &orig,
                     Eigen::MatrixXd &Vc, bool remove_interior_vertices) {
  Eigen::VectorXi new_orig = orig;

  if (remove_interior_vertices) {
    new_orig = remove_interior(F, orig);
  }

  // Get the correct indices.
  Eigen::VectorXi b = extract_fixed_vertices(orig);

  // The positions of these indices are held constant (just keep the input values)
  Eigen::MatrixXd V_bc = igl::slice(V, b, 1);

  Eigen::SparseMatrix<double> L, M;
  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT,M);

  // Initialize Vc with input vertices V
  Vc = V;

  igl::viewer::Viewer v;
  int harmonic_idx = 2;
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
    } else if (key == '1') {
      printf("Setting harmonic to be 1\n");
      harmonic_idx = 1;
    } else if (key == '1') {
      printf("Setting harmonic to be 2\n");
      harmonic_idx = 2;
    } else if (key == '3') {
      printf("Setting harmonic to be 3\n");
      harmonic_idx = 3;
    } else if (key == ' ') {
      // Press spacebar to get the next version.
      
      // Recompute M, leave L alone.
      igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
      //igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_VORONOI,M);
      //igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);

      // How much should we change by?
      Eigen::MatrixXd D;
      //Eigen::MatrixXd V_bc = igl::slice(Vc, b, 1);
      bool success = igl::harmonic(L,M, b,V_bc, harmonic_idx, D);
      printf("SUCCESS: %d\n", success);
      //igl::harmonic(Vc,F, b,V_bc, 2, D);

      cout << "lengths" << endl;
      cout << Vc.rows() << " " << D.rows() << endl;
      // Calculate the difference.
      double diff = 0;
      // Vertices updated like:
      for (int i = 0; i < V.rows(); ++i) {
        double tmp = (Vc.row(i) - D.row(i)).norm();
        // if (std::isnan(tmp)) {
        //   cout << "index: " << i << endl;
        //   cout << "Vc: " << endl;
        //   cout << Vc.row(i) << endl;
        //   cout << "D: " << endl;
        //   cout << D.row(i) << endl;
        //   cout << endl;
        // }

        diff += (Vc.row(i) - D.row(i)).norm();
      }
      printf("Difference this round is %lf\n", diff);
      
      // Update the values.
      //Vc = Vc + D;
      Vc = D;

      // Calculate the energy.
      Eigen::SparseMatrix<double> M_new, Mi, L_new;
      igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_DEFAULT,M_new);
      igl::invert_diag(M_new,Mi);
      igl::cotmatrix(Vc,F,L_new);
      auto en = (Vc.transpose() * L_new * Mi * L_new * Vc).trace();
      // Energy needs to be computed on the new Laplacian matrix
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
  printf("  Press 1,2,or 3 to use different versions of the harmonic.\n");
  printf("  Close down image to return shape to previous function\n");
  v.data.set_mesh(Vc, F);
  Eigen::MatrixXd cols;
  igl::jet(orig, true, cols);
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
