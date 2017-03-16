#include "curvatureFlow.h"

#include <iostream> // cout
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
#include <igl/remove_unreferenced.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "Helpers.h"

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

// Will add the faces to border_F
void getNewPathPoints(const Eigen::MatrixXd &V,
                      const std::vector<std::vector<int> > &adj_list,
                      std::map<int, std::map<int, int> > &edge_count,
                      std::set<int> &visited,
                      int start, int offset, bool pos,
                      std::vector<Eigen::VectorXi> &border_F) {
  // Create an ordering of top_v around the mesh.
  std::vector<int> path;

  // Start with one top_v_orig
  path.push_back(start);
  visited.insert(start);

  int next = -2;

  while(next != -1) {
    // Add them to the path
    int current = (next == -2) ? start : next;
    next = -1;

    for (int n : adj_list[current]) {
      // Haven't found this point, so go there.
      
      if (visited.find(n) == visited.end() && edge_count[current][n] == 1) {
        next = n;
        path.push_back(next);
        visited.insert(next);
        break;
      }
    }
  }

  std::vector<Eigen::VectorXi> temp_f;
  // Add two new faces for each vertex in the path.
  for (int i = 0; i < path.size(); ++i) {
    int next = (i + 1) % path.size();

    // Create two new faces.
    Eigen::Vector3i f1;
    f1 << offset + path[i], offset + path[next], path[i];
    Eigen::Vector3i f2;
    f2 << path[next], path[i], offset + path[next];
    temp_f.push_back(f1);
    temp_f.push_back(f2);
  }

  const auto &p1 = V.row(path[0]);
  const auto &p2 = V.row(path[1]);
  const auto &p3 = V.row(path[2]);
  Eigen::Vector3d u = (p3 - p2).normalized();
  Eigen::Vector3d v = (p1 - p2).normalized();
  Eigen::Vector3d n = u.cross(v);
  //printf("n2:%lf pos:%d\n", n(2), pos);
  int numGreater = 0, numLess = 0;
  for (int i = 0; i < path.size(); ++i) {
    const auto &a = V.row(path[i]);
    const auto &b = V.row(path[ (i+1) % path.size()]);
    const auto &c = V.row(path[ (i+2) % path.size()]);
    Eigen::Vector3d u2 = (c - b);
    Eigen::Vector3d v2 = (a - b);
    Eigen::Vector3d n2 = u2.cross(v2);
    if (n2(2) > GLOBAL::EPS) {
      numGreater++;
    } else if (n2(2) < GLOBAL::EPS) {
      numGreater--;
    }
  }
  if ((numGreater > numLess) == pos) {
    // Flip face normals.
    for (auto &f : temp_f) {
      int temp = f(0);
      f(0) = f(1);
      f(1) = temp;
    }
  }

  // Add them all to border_F
  for (auto &f : temp_f) {
    border_F.push_back(f);
  }
}

void duplicateTopVertices(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &orig,
                          Eigen::MatrixXd &newV, Eigen::MatrixXi &newF,
                          Eigen::VectorXi &newOrig, 
                          bool use_max,
                          double ext_amt=0.1) {
  // Identify extreme vertices
  std::vector<int> top_v_orig;
  std::set<int> top_v_nonorig;
  std::set<int> top_v_all;
  Eigen::VectorXd limits;
  if (use_max) {
    limits = V.colwise().maxCoeff();
  } else {
    limits = V.colwise().minCoeff();
  }

  for (int i = 0; i < V.rows(); ++i) {
    if (std::abs(V(i, 2) - limits(2)) < GLOBAL::EPS) {
      if (orig(i) != GLOBAL::nonoriginal_marker) {
        top_v_orig.push_back(i);
      } else {
        top_v_nonorig.insert(i);
      }

      top_v_all.insert(i);
    }
  }

  // Degenerate case when there is only one maximal point.
  if (top_v_orig.size() < 1) {
    newV = V;
    newF = F;
    newOrig = orig;
    return;
  }

  // edge_count[u][v] = edge_count[v][u]
  std::map<int, std::map<int, int> > edge_count;
  for (int i = 0; i < F.rows(); i++) {
    bool do_it = true;
    for (int j = 0; j < 3; j++) {
      int u = F(i, j);
      if (top_v_all.find(u) == top_v_all.end())
        do_it = false;
    }

    if (!do_it)
      continue;

    for (int j = 0; j < 3; j++) {
      int u = F(i, j);
      int v = F(i, (j+1) % 3);

      edge_count[u][v]++;
      edge_count[v][u]++;
    }
  }

  int offset = V.rows();
  // Duplicate points (kinda overkill)
  Eigen::MatrixXd add_V(V.rows() * 2, V.cols());
  for (int i = 0; i < V.rows(); ++i) {
    add_V.row(i) = V.row(i);
    add_V.row(offset + i) = V.row(i);
    if (use_max) {
      add_V(offset + i, 2) += ext_amt;
    } else {
      add_V(offset + i, 2) -= ext_amt;
    }
  }


  // First, create an adjacency list.
  std::vector<std::vector<int> > adj_list;
  igl::adjacency_list(F, adj_list);
  std::set<int> visited; // vertices we've visited.
  // Extra faces we should add.
  std::vector<Eigen::VectorXi> border_F;
  // Use each vertex as a starting point; won't add anything if we've already
  // visited.
  for (int idx : top_v_orig) {
    if (visited.find(idx) == visited.end()) {
      getNewPathPoints(V, adj_list, edge_count, visited, idx, offset, use_max, border_F);
    }
  }

  Eigen::MatrixXi add_F(F.rows() + border_F.size(), F.cols());

  int numFound = 0;
  // Look at each of the non-shell points to see if we should raise this face.
  for (int i = 0; i < F.rows(); ++i) {
    add_F.row(i) = F.row(i);

    bool found = false;
    for (int j = 0; j < F.cols(); ++j) {
      // New faces from new points.
      if (top_v_nonorig.find(F(i, j)) != top_v_nonorig.end()) {
        found = true;
        break;
      }
    }
    // Also allowed if all 3 faces are in top_v_all (all vertices on the top)
    if (!found) {
      if (top_v_all.find(F(i, 0)) != top_v_all.end() &&
          top_v_all.find(F(i, 1)) != top_v_all.end() &&
          top_v_all.find(F(i, 2)) != top_v_all.end()) {
        found = true;
      }
    }

    // If we've found the vertex, add new faces.
    if (found) {
      numFound++;
      for (int j = 0; j < F.cols(); ++j) {
        add_F(i, j) += offset; // it has a new vertex.
      }
    }
  }
  /*
  printf("Found %d faces to change, size of nonorig is %d, orig is %d\n", 
         numFound, top_v_nonorig.size(), top_v_orig.size());
         */

  // Add the additional faces.
  for (int i = 0; i < border_F.size(); ++i) {
    add_F.row(F.rows() + i) = border_F[i];
  }

  // Now, remove unreferenced vertices.
  Eigen::VectorXi I, J; // J is new indices into F, I new indices into V
  igl::remove_unreferenced(add_V.rows(), add_F, J, I);
  //printf("J %lu, I %lu, V: %lu F: %lu\n", J.rows(), I.rows(), add_V.rows(), add_F.rows());

  // Update F
  newF = add_F;
  for (int i = 0; i < newF.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      // Get the new vertex.
      newF(i, j) = J(newF(i, j));
    }
  }
  // Update V, with only things from I.
  igl::slice(add_V, I, 1, newV);

  newOrig.resize(newV.rows());
  for (int i = 0; i < newOrig.rows(); ++i) {
    int idx = I(i);
    // If it's an added point
    if (idx >= V.rows()) {
      // Get the original offset marker.
      newOrig(i) = orig(idx - offset);
    } else {
      newOrig(i) = orig(idx);
    }
  }
}

} // namespace

double biharmonic_new(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    const Eigen::VectorXi &orig,
                    Eigen::MatrixXd &Vc, Eigen::MatrixXi &Fc,
                    Eigen::VectorXi &Mc,
                    double change_val) {
  Eigen::VectorXi new_orig = orig;
  removeInterior(F, orig, new_orig);
  Eigen::MatrixXd newV, V2;
  Eigen::MatrixXi F2;
  Eigen::VectorXi M2;
  //printf("inside biharmonic_new\n");
  //Helpers::viewTriMesh(V, F, new_orig);
  duplicateTopVertices(V, F, new_orig, V2, F2, M2, true); // above
  //printf("After extending top\n");
  //Helpers::viewTriMesh(V2, F2, M2);
  duplicateTopVertices(V2, F2, M2, newV, Fc, Mc, false); // below.
  //printf("Before biharmonic\n");
  //Helpers::viewTriMesh(newV, Fc, Mc);
  return biharmonic(newV, Fc, Mc, Vc, false /* don't remove interior */);
}

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
  Eigen::SparseMatrix<double> Lt;
  igl::cotmatrix(Vc,F,Lt);
  igl::massmatrix(Vc,F,igl::MASSMATRIX_TYPE_DEFAULT,M);
  Eigen::SparseMatrix<double> Mi;
  igl::invert_diag(M, Mi);
  return (Vc.transpose() * Lt * Mi * Lt * Vc).trace();
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
  //convertVectorToIndices(new_orig, b);
  convertVectorToIndices(orig, b);
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
      igl::harmonic(L,M, b,V_bc, harmonic_idx, D);
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
