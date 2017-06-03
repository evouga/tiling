#include "offsetSurface.h"

#include <map>
#include <utility>
#include <vector>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <igl/adjacency_list.h>
#include <igl/boundary_facets.h>
#include <igl/collapse_edge.h>
#include <igl/copyleft/cgal/complex_to_mesh.h>
#include <igl/edge_flaps.h>
#include <igl/exterior_edges.h>
#include <igl/in_element.h>
#include <igl/jet.h>
#include <igl/remove_unreferenced.h>
#include <igl/unique.h>
#include <igl/unique_edge_map.h>
#include <igl/viewer/Viewer.h>

#include "glob_defs.h"
#include "Helpers.h"
#include "marching_tets.h"

#define OFFSET_DEBUG 0

using namespace std;

// begin namespace for helper methods.
namespace {

bool should_collapse(int a, int b, double zmin, double zmax,
                     const Eigen::MatrixXd &V,
                     const Eigen::VectorXi &M) {
  // If a is original and b is not and both are on the boundary, return true
  if (M(a) != GLOBAL::nonoriginal_marker &&
      M(b) == GLOBAL::nonoriginal_marker) {
    return (std::fabs(V(a, 2) - zmin) < GLOBAL::EPS &&
            std::fabs(V(b, 2) - zmin) < GLOBAL::EPS) ||
           (std::fabs(V(a, 2) - zmax) < GLOBAL::EPS &&
            std::fabs(V(b, 2) - zmax) < GLOBAL::EPS);
  }
  // Otherwise, return false.
  return false;
}

bool should_collapse_type(int a, int b, double zmin, double zmax,
                          const Eigen::MatrixXd &V,
                          const Eigen::VectorXi &M,
                          const Eigen::VectorXi &VValence,
                          int type) {
  if (!should_collapse(a,b, zmin,zmax, V,M)) return false;

  // Know that one of them must be nonoriginal
  if (type == 1) {
    if (M(a) != GLOBAL::nonoriginal_marker) {
      return VValence(a) == 3;
    }
    return VValence(b) == 3;
  } 
  // Type 2 is return true if outer ring has valence 3.
  if (type == 2) {
    if (M(a) == GLOBAL::nonoriginal_marker) {
      return VValence(a) == 3;
    }
    return VValence(b) == 3;
  }
  // Type 0 is all others.
  return true;
}

void getVertexValence(const Eigen::MatrixXd &V,
                      const Eigen::MatrixXi &F,
                      Eigen::VectorXi &VValence) {
  VValence.resize(V.rows());
  VValence.setZero();
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < F.cols(); ++j) {
      // Ignore faces that have already been collapsed.
      if (F(i, j) == IGL_COLLAPSE_EDGE_NULL) break;

      VValence(F(i, j))++;
    }
  }
}

vector<int> getEdgesToRemove(const Eigen::MatrixXd &V,
                             const Eigen::MatrixXi &E,
                             const Eigen::VectorXi &M) {
  // The boundary vertices.
  double z_min = 1e10;
  double z_max = -1e10;

  for (int i = 0; i < V.rows(); i++) {
    z_min = min(z_min, V(i, 2));
    z_max = max(z_max, V(i, 2));
  }

  vector<int> to_remove;
  // Check all the edges for things we should collapse.
  for (int i = 0; i < E.rows(); ++i) {
    // Check for both directions.
    if (should_collapse(E(i, 0), E(i, 1), z_min, z_max, V, M) ||
        should_collapse(E(i, 1), E(i, 0), z_min, z_max, V, M)) {
      if (E(i, 0) == E(i, 1)) {
        printf("ERROR in computing!!!\n");
      }
      if (E(i, 0) == 21 || E(i, 0) == 19 || i == 100) {
        printf("Here, we have %d: %d=%d and %d=%d\n", 
               i, E(i, 0), M(E(i, 0)), E(i, 1), M(E(i, 1)));
      }
      to_remove.push_back(i);
    }
  }

  return to_remove;
}

// Returns list of (a, b) pairs, where "a" is the vertex to keep.
vector<pair<int, int> > getEdgesToRemove(const Eigen::MatrixXd &V,
                                         const vector<vector<int> > &graph,
                                         const Eigen::MatrixXi &O) {
  // The boundary vertices.
  double z_min = 1e10;
  double z_max = -1e10;

  for (int i = 0; i < V.rows(); i++) {
    z_min = min(z_min, V(i, 2));
    z_max = max(z_max, V(i, 2));
  }

  // Create [a, b] pairs, where a is on the contour.
  vector<pair<int, int> > to_remove;

  for (int a = 0; a < graph.size(); a++) {
    // If a is not on the contour, ignore it.
    if (O(a) == GLOBAL::nonoriginal_marker)
      continue;

    for (int b : graph[a]) {
      // If b is on the contour, ignore it.
      if (O(b) != GLOBAL::nonoriginal_marker)
        continue;

      double z = V(b, 2);

      // Vertex b lies on the boundary.
      if (abs(z - z_min) <= 1e-5 || abs(z - z_max) <= 1e-5)
        to_remove.push_back(make_pair(a, b));
    }
  }

  return to_remove;
}

// Collapse vertex b onto vertex a.
bool removeEdge(int a, int b, Eigen::MatrixXi &F,
                vector<vector<int> > &vertex_to_faces) {
  bool found_b = false;

  // Replace all occurrences of b with a.
  for (int face_index : vertex_to_faces[b]) {
    for (int i = 0; i < F.cols(); i++) {
      if (F(face_index, i) == b) {
        F(face_index, i) = a;
        found_b = true;
      }
    }
  }

  // No vertex b's were found, vertex was already removed.
  if (!found_b)
    return false;

  int number_triangles_collapsed = 0;

  // Now, the two triangles straddling the edge a-b will have [a, a, x].
  for (int face_index : vertex_to_faces[b]) {
    int number_a = 0;

    for (int i = 0; i < F.cols(); i++) {
      if (F(face_index, i) == a)
        number_a++;
    }

    // Set straddle triangles to [-1, -1, -1].
    if (number_a == 2) {
      number_triangles_collapsed++;

      for (int i = 0; i < F.cols(); i++)
        F(face_index, i) = -1;
    }
  }

  return (number_triangles_collapsed == 2);
}

// A collapsed face is [-1, -1, -1].
bool isDeadFace(const Eigen::RowVectorXi &face) {
  for (int i = 0; i < 3; i++) {
    if (face(i) != -1)
      return false;
  }
  return true;
}

// Removes all faces that are [-1, -1, -1].
void cleanupDeadFaces(const Eigen::MatrixXd &V_old, const Eigen::MatrixXi &F_old,
                      const Eigen::VectorXi &O_old,
                      Eigen::MatrixXd &V_new, Eigen::MatrixXi &F_new,
                      Eigen::VectorXi &O_new) {
  // Find number alive.
  int number_alive = 0;
  for (int i = 0; i < F_old.rows(); i++) {
    if (isDeadFace(F_old.row(i)))
      continue;

    number_alive++;
  }

  // Remove the dead faces.
  Eigen::MatrixXi F_prepared(number_alive, 3);
  int next_spot = 0;
  for (int i = 0; i < F_old.rows(); i++) {
    if (isDeadFace(F_old.row(i)))
      continue;

    F_prepared.row(next_spot++) = F_old.row(i);
  }

  // Mapping from old index to new index. New index is -1 if removed.
  Eigen::VectorXi I;
  igl::remove_unreferenced(V_old, F_prepared, V_new, F_new, I);

  // Update markers.
  O_new.resize(V_new.rows());
  O_new.setZero();

  for (int i = 0; i < O_old.rows(); i++) {
    // Vertex was removed.
    if (I(i) == -1)
      continue;

    O_new(I(i)) = max(O_new(I(i)), O_old(i));
  }
}

// Remove N1 neighbors of contour vertices.
//
// type is one of:
//   0: collapse any edge
//   1: collapse ring vertices with valence three
//   2: collapse outer vertices with valence three
int removeRing(const Eigen::MatrixXd &V_old, const Eigen::MatrixXi &F_old,
               const Eigen::VectorXi &O_old,
               Eigen::MatrixXd &V_new, Eigen::MatrixXi &F_new,
               Eigen::VectorXi &O_new, int type) {
  // Start of new code.
  Eigen::MatrixXi E, EF,EI;
  Eigen::VectorXi EMAP;
  igl::edge_flaps(F_old, E,EMAP,EF,EI);

  // Also get vertex valences
  Eigen::VectorXi VValence;
  getVertexValence(V_old, F_old, VValence);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi O;

  V = V_old;
  F = F_old;
  O = O_old;
  // The boundary vertices.
  double z_min = 1e10;
  double z_max = -1e10;

  for (int i = 0; i < V.rows(); i++) {
    z_min = min(z_min, V(i, 2));
    z_max = max(z_max, V(i, 2));
  }

  int num_removed = 0;
  std::set<int> collapsed;
  /*
  igl::viewer::Viewer viewer;
  Eigen::MatrixXd Cols;
  igl::jet(O, true, Cols);
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(Cols);
  int i = 0;
  viewer.callback_key_down = [&](igl::viewer::Viewer& vv,
                                 unsigned char key, int modifier) {
    if (key == ' ') {
      if (i >= E.rows()) {
        printf("Can't proceed, try quitting.\n");
        return false;
      }
      int prev_idx_a = -1, prev_idx_b = -1;
      int prev_val_a, prev_val_b;
      int new_val;
      for (;i < E.rows(); ++i) {
        // Means it's already done.
        if (E(i, 0) == E(i, 1)) {
          //printf("Already done edge with %d=%d\n", E(e, 0), E(e, 1));
          continue;
        }
        // Only use each vertex once.
        if (collapsed.find(E(i, 0)) != collapsed.end() ||
            collapsed.find(E(i, 1)) != collapsed.end())
          ;//continue;
        // Ignore if it's not valid.
        if (!should_collapse(E(i, 0), E(i, 1), z_min, z_max, V, O) &&
            !should_collapse(E(i, 1), E(i, 0), z_min, z_max, V, O)) {
          continue;
        }

        int e = i;
        // Get the index of the original one.
        int orig = E(e, 0), nonorig = E(e, 1);
        if (O(orig) == GLOBAL::nonoriginal_marker) {
          orig = E(e, 1);
          nonorig = E(e, 0);
        }
        if (O(orig) == GLOBAL::nonoriginal_marker) {
          printf("[%d] This is a super big error, with %d=%d, %d=%d\n", 
                 e, orig,O(orig), nonorig,O(nonorig));
        }
        bool worked = igl::collapse_edge(e, V.row(orig), V,F,E, EMAP, EF,EI);
        if (!worked) {
          printf("Error: edge collapse didn't work for edge %d=%d,%d (orig:%d,%d)\n",
                 e, E(e, 0), E(e, 1), orig, nonorig);
          //O(orig) = e;
          //O(nonorig) = e;
          new_val = 10;
        } else {
          printf("       Edge collapse DID WORK!!! for edge %d=%d,%d (orig:%d,%d)\n",
                 e, E(e, 0), E(e, 1), orig, nonorig);
          new_val = 9;
          // Also need to update the original value.
          O(nonorig) = O(orig);
          collapsed.insert(orig);
          collapsed.insert(nonorig);
          num_removed++;
        }
        prev_idx_a = orig;
        prev_idx_b = nonorig;
        i++;
        break;
      }

      if (prev_idx_a != -1) {
        prev_val_a = O(prev_idx_a);
        prev_val_b = O(prev_idx_b);
        O(prev_idx_a) = new_val;
        O(prev_idx_b) = new_val;
        igl::jet(O, true, Cols);
        vv.data.clear();
        vv.data.set_mesh(V, F);
        vv.data.set_colors(Cols);
        O(prev_idx_a) = prev_val_a;
        O(prev_idx_b) = prev_val_b;
        return true;
      }
      return false;
    }
  };
  viewer.launch();
  */

  // Get edges to remove.
  //vector<int> remove = getEdgesToRemove(V, E, O);
  //for (int e : remove) {
  for (int i = 0; i < E.rows(); ++i) {
    printf("i is %d/%d (%d,%d)\n", i, E.rows(), EF(i, 0), EF(i, 1));
    // Means it's already done.
    if (E(i, 0) == E(i, 1)) {
      //printf("Already done edge with %d=%d\n", E(e, 0), E(e, 1));
      continue;
    }
    // Only use each vertex once.
    if (collapsed.find(E(i, 0)) != collapsed.end() ||
        collapsed.find(E(i, 1)) != collapsed.end())
      continue;

    if (EF(i, 0) < 0 || EF(i, 1) < 0) {
      continue;
    }
    // Face check, already collapsed some face.
    if (F(EF(i, 0), 0) == IGL_COLLAPSE_EDGE_NULL ||
        F(EF(i, 1), 0) == IGL_COLLAPSE_EDGE_NULL) {
      printf("Found bad face!!!\n");
      continue;
    }

    // Ignore if it's not valid.
    if (!should_collapse_type(E(i, 0), E(i, 1), z_min, z_max, V, O, VValence, type) &&
        !should_collapse_type(E(i, 1), E(i, 0), z_min, z_max, V, O, VValence, type)) {
      continue;
    }

    int e = i;
    printf("e is %d(%d,%d)\n", e, E(e, 0), E(e, 1));
    if (e == 107) {
      igl::viewer::Viewer viewer;
      Eigen::VectorXi Ot = O;
      Ot(E(e, 0)) = 10;
      Ot(E(e, 1)) = 10;
      Eigen::MatrixXd Cols;
      igl::jet(Ot, true, Cols);
      viewer.data.set_mesh(V, F);
      viewer.data.set_colors(Cols);
      viewer.launch();
    }

    // Get the index of the original one.
    int orig = E(e, 0), nonorig = E(e, 1);
    if (O(orig) == GLOBAL::nonoriginal_marker) {
      orig = E(e, 1);
      nonorig = E(e, 0);
    }
    if (O(orig) == GLOBAL::nonoriginal_marker) {
      printf("[%d] This is a super big error, with %d=%d, %d=%d\n", 
             e, orig,O(orig), nonorig,O(nonorig));
    }
    printf("collapsing %d to (%lf,%lf,%lf)... (%d,%d,%d) (%d,%d,%d)\n", e, 
           V(orig, 0), V(orig, 1), V(orig, 2),
           F(EF(e, 0), 0), F(EF(e, 0), 1), F(EF(e, 0), 2),
           F(EF(e, 1), 0), F(EF(e, 1), 1), F(EF(e, 1), 2));
    fflush(stdout);
    bool worked = igl::collapse_edge(e, V.row(orig), V,F,E, EMAP, EF,EI);
    printf("done!\n");
    fflush(stdout);
    if (!worked) {
      printf("Error: edge collapse didn't work for edge %d=%d,%d (orig:%d,%d)\n",
             e, E(e, 0), E(e, 1), orig, nonorig);
      //O(orig) = e;
      //O(nonorig) = e;
    } else {
      printf("       Edge collapse DID WORK!!! for edge %d=%d,%d (orig:%d,%d)\n",
             e, E(e, 0), E(e, 1), orig, nonorig);
      num_removed++;
      // Also need to update the original value.
      O(nonorig) = O(orig);
      collapsed.insert(orig);
      collapsed.insert(nonorig);
    }
  }
  printf("Finished collapsing everything...\n");

  // Get valid faces.
  vector<int> valid_f;
  for (int i = 0; i < F.rows(); ++i) {
    int valid = 0;
    for (int j = 0; j < F.cols(); ++j) {
      valid++;
      if (F(i, j) == IGL_COLLAPSE_EDGE_NULL) {
        valid--;
      }
    }
    if (valid == 3) {
      valid_f.push_back(i);
    } else if (valid > 0) {
      printf("Found a face with %d valid vertices?\n", valid);
    }
  }
  
  F_new.resize(valid_f.size(), 3);
  for (int i = 0; i < valid_f.size(); ++i) {
    F_new.row(i) = F.row(valid_f[i]);
  }
  V_new = V;
  O_new = O;
    
  return num_removed;
  //
  // End of new code.
  //
  //
  //
  //

  vector<vector<int> > graph;
  igl::adjacency_list(F_old, graph);

  // List of (keep, remove) pairs of vertices.
  vector<pair<int, int> > to_remove = getEdgesToRemove(V_old, graph, O_old);


  // Mapping of vertex index to list of face indices it occurs in.
  vector<vector<int> > vertex_to_faces(V_old.rows());
  for (int i = 0; i < F_old.rows(); i++) {
    for (int j = 0; j < F_old.cols(); j++)
      vertex_to_faces[F_old(i, j)].push_back(i);
  }

  Eigen::MatrixXi F_with_dead_faces = F_old;

  // Decimate the edges.
  set<int> have_removed;
  for (auto a_b : to_remove) {
    // Want to get rid of vertex b.
    int a = a_b.first;
    int b = a_b.second;

    // Have already removed b.
    if (have_removed.find(b) != have_removed.end())
      continue;

    // Collapsed faces will be [-1, -1, -1].
    if (removeEdge(a, b, F_with_dead_faces, vertex_to_faces)) {
      //printf("Collapsed edge %d-%d\n", b, a);
      have_removed.insert(b);
    }
  }

  cleanupDeadFaces(V_old, F_with_dead_faces, O_old,
                   V_new, F_new, O_new);

  return (V_old.rows() - V_new.rows());
}
  
} // end namespace for helper methods.

namespace OffsetSurface {

// Needed CGAL typedefs
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 CTr;
typedef CGAL::Complex_2_in_triangulation_3<CTr> C2t3;
typedef CTr::Geom_traits          GT;
typedef GT::Sphere_3              Sphere_3;
typedef GT::Point_3               Point_3;
typedef GT::FT                    FT;
//typedef FT (*Function)(Point_3);
typedef std::function<FT (Point_3)> Function;
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

/*
 * Function to call when you have a triangulation.
 */
void generateOffsetSurface(Triangulation &T, double off,
                           Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff) {
  const igl::AABB<Eigen::MatrixXd,3> &tree = T._tree;
  const Eigen::MatrixXd &V = T._V;
  const Eigen::VectorXd &C = T._C;
  const Eigen::MatrixXi &TT = T._T; // tets
  auto maxs = V.colwise().maxCoeff();
  auto mins = V.colwise().minCoeff();

  // Need the triangulation and the offset
  Function fun = 
      [&](const Point_3 &qq) -> FT
      {
        Eigen::RowVectorXd q(3);
        q << qq.x(), qq.y(), qq.z();
        // Only return first element containing q
        std::vector<int> tets = tree.find(V, TT, q, true);
        // Make sure it's in a tet
        if (tets.size() < 1) {
          return GLOBAL::highest_temp;
        }
        int tet = tets[0];

        // Get the 4 vertices
        Eigen::MatrixXd mat(4, 4);
        bool z_boundary_bad = false;
        for (int i = 0; i < 4; ++i) {
          mat(0, i) = V(TT(tet,i), 0);
          mat(1, i) = V(TT(tet,i), 1);
          mat(2, i) = V(TT(tet,i), 2);
          mat(3, i) = 1;
          if (OFFSET_DEBUG) 
            printf("%f %f %f %f 1\n",
                   V(TT(tet,i), 0), V(TT(tet,i), 1), V(TT(tet,i), 2),
                   C(TT(tet,i)));
          // If one of the vertices is on the z boundary, and this is also on
          // the z boundary, then it's not a real point.
          if ( (V(TT(tet, i), 2) == maxs(2) ||
                V(TT(tet, i), 2) == mins(2)) &&
               C(TT(tet, i)) == GLOBAL::outside_temp) {
            z_boundary_bad = true;
          }
        }
        // If the z-coordinate is on the boundary and one of the tet
        // vertices is 'outside', return a bad value as well.
        if (z_boundary_bad && (q(2) == maxs(2) ||
                               q(2) == mins(2))) {
          return GLOBAL::highest_temp;
        }

        Eigen::Vector4d query(4);
        query(0) = q(0);
        query(1) = q(1);
        query(2) = q(2);
        query(3) = 1;

        Eigen::Vector4d alpha = mat.inverse() * query;
        double val = alpha(0)*C(TT(tet, 0)) +
                     alpha(1)*C(TT(tet, 1)) +
                     alpha(2)*C(TT(tet, 2)) +
                     alpha(3)*C(TT(tet, 3));
        /*
        std::cout << "Mult is:\n"
             << mat << "\n"
             << mat.inverse() << "\n"
             << query << "\n"
             << alpha << "\n"
             << val << "\n"
             << C(TT(tet, 0)) << ", " << C(TT(tet, 1)) << ", " 
             << C(TT(tet, 2)) << ", " << C(TT(tet, 3)) << "\n";
             */
        double col_val = 0;
        for (int i = 0; i < 4; ++i) {
          col_val += V(TT(tet, i), 0) +
              V(TT(tet, i), 1) +
              V(TT(tet, i), 2);
        }
        if (OFFSET_DEBUG) printf("%f %f %f %f 0\n",
                                 q(0), q(1), q(2), val);
        if (OFFSET_DEBUG) printf("%f %f %f %f 2\n",
                                 q(0), q(1), q(2), col_val);
        return val - off;
      };

  /*
  for (double x = -0.51; x <= 0.51; x += 0.02) {
    for (double y = -0.51; y <= 0.51; y += 0.02) {
      for (double z = -0.51; z <= 0.51; z += 0.02) {
        Point_3 q(x, y, z);
        fun(q);
      }
    }
  }
  fprintf(stderr, "Finished with grid!\n");
  return;
  */
  // Code borrowed from igl: https://github.com/libigl/libigl/blob/master/include/igl/copyleft/cgal/signed_distance_isosurface.cpp
  const auto &Vmax = V.colwise().maxCoeff();
  const auto &Vmin = V.colwise().minCoeff();
  const double bbd = (Vmax-Vmin).norm();
  const double r = bbd/2.;
  const auto &Vmid = 0.5*(Vmax + Vmin);

  Point_3 cmid(Vmid(0), Vmid(1), Vmid(2));
  //Point_3 cmid(0, 0, 0);
  // Now, find the surface
  //Sphere_3 bounding_sphere(cmid, r*r);  // squared radius
  Sphere_3 bounding_sphere(cmid, (r+off)*(r+off));  // squared radius
  Surface_3 surface(fun, bounding_sphere);
  // Use default values for angle bound, radius bound, distance bound (respectively)
  //CGAL::Surface_mesh_default_criteria_3<CTr> criteria(28, 0.2, 0.5);
  //CGAL::Surface_mesh_default_criteria_3<CTr> criteria(28, 1.1, 1.1);
  CGAL::Surface_mesh_default_criteria_3<CTr> criteria(28, 2.1, 2.1);
  // Mesh surface
  printf("Inserting points to triangulation...");
  CTr tr;  // 3D-Delaunay triangulation
  for (int i = 0; i < V.rows(); ++i) {
    tr.insert(Point_3(V(i,0), V(i,1), V(i,2)));
  }
  printf("done!\nCreating surface mesh...");
  fflush(stdout);
  C2t3 c2t3(tr); // 2D-complex in 3D-Delaunay triangulation
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  printf("done!\nNow converting complex to mesh...");
  fflush(stdout);
  // Complex to (V,F)
  igl::copyleft::cgal::complex_to_mesh(c2t3, Voff, Foff);
  printf("done!\n");
}

/**
 * Initial function call, 'returns' the triangulation for successive calls
 */
void generateOffsetSurface(const Eigen::MatrixXd &V,
                           const Eigen::MatrixXi &TT, const Eigen::VectorXd &C,
                           double off, Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff,
                           Triangulation &T) {
  igl::AABB<Eigen::MatrixXd,3> &tree = T._tree;
  tree.init(V, TT);
  T._V = V;
  T._T = TT;
  T._C = C;

  generateOffsetSurface(T, off, Voff, Foff);
}

bool fitsOffset(int idx, double off, const Eigen::MatrixXi &TT, const Eigen::VectorXd &C) {
  for (int i = 0; i < TT.cols(); ++i) {
    if (C(TT(idx, i)) > off) {
      return false;
    }
  }
  return true;
}
bool hasOriginal(int idx, const Eigen::MatrixXi &TT, const Eigen::VectorXi &TO) {
  for (int i = 0; i < TT.cols(); ++i) {
    if (TO(TT(idx, i)) != GLOBAL::nonoriginal_marker) {
      return true;
    }
  }
  return false;
}
bool hasBoundary(int idx, const Eigen::MatrixXi &TT, const Eigen::VectorXd &C) {
  for (int i = 0; i < TT.cols(); ++i) {
    if (C(TT(idx, i)) == GLOBAL::outside_temp) {
      return true;
    }
  }
  return false;
}
void generateOffsetSurface_naive(const Eigen::MatrixXd &V,
                                 const Eigen::MatrixXi &TT,
                                 const Eigen::VectorXi &TO,
                                 const Eigen::VectorXd &C, double off,
                                 Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff,
                                 Eigen::VectorXi &Ooff) {
  std::vector<int> s;
  for (int i = 0; i < TT.rows(); ++i) {
    // Either all tet corners are within the offset
    // or there is an original vertex and no boundary (outside_temp) vertices
    if (fitsOffset(i, off, TT, C) ||
        (hasOriginal(i, TT, TO) && !hasBoundary(i, TT, C))) {
      s.push_back(i);
    }
  }

  Voff.resize(s.size() * 4, 3);
  Foff.resize(s.size() * 4, 3);
  Ooff.resize(s.size() * 4);

  for (unsigned i = 0; i < s.size(); ++i)
  {
    for (int j = 0; j < 4; ++j) {
      Voff.row(i*4+j) = V.row(TT(s[i],j));
    }
    Foff.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
    Foff.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
    Foff.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
    Foff.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
    for (int j = 0; j < 4; ++j) {
      Ooff(i*4+j) = TO(TT(s[i],j));
    }
  }

  // Make sure the vertices are all unique
	// Get all the unique vertices
	Eigen::MatrixXd Vu;
  Eigen::VectorXi Ou;
	// contains mapping from unique to all indices
	// Size: #V
	Eigen::VectorXi unique_to_all;
	// contains mapping from all indices to unique
	// Size: #V_rep
	Eigen::VectorXi all_to_unique;
	igl::unique_rows(Voff, Vu,unique_to_all,all_to_unique);
  // Remember these.
  Voff = Vu;

  // Update the original markers.
  Ou.resize(Voff.rows());
  Ou.setZero();
  for (int i = 0; i < Ou.rows(); ++i) {
    Ou(i) = Ooff(unique_to_all(i));
  }
  Ooff = Ou;

  // Also need to update faces.
  for (int i = 0; i < Foff.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      Foff(i, j) = all_to_unique(Foff(i, j));
    }
  }
  // And remove duplicate faces (creates a shell)
  Eigen::MatrixXi Fu;
  igl::unique_rows(Foff, Fu,unique_to_all,all_to_unique);
  Foff = Fu;
}

// Will push the nonoriginal vertices on the boundary up or down.
void quickHack(Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXi &M) {
  // The boundary vertices.
  double zmin = 1e10;
  double zmax = -1e10;
  for (int i = 0; i < V.rows(); i++) {
    zmin = min(zmin, V(i, 2));
    zmax = max(zmax, V(i, 2));
  }

  // Adjust everything up or down.
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < F.cols(); ++j) {
      // If it's not on the boundary, don't care.
      if (V(F(i, j), 2) != zmin && V(F(i, j), 2) != zmax) continue;
      // Don't touch things that are original.
      if (M(F(i, j)) != GLOBAL::nonoriginal_marker) continue;
      // Move everything up by epsilon.
      if (V(F(i, j), 2) == zmin) {
        V(F(i, j), 2) += 2 * GLOBAL::EPS;
        M(F(i, j)) = 10;
      } else {
        V(F(i, j), 2) -= 2 * GLOBAL::EPS;
        M(F(i, j)) = 9;
      }

    }

  }
}

void marchingOffsetSurface(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO,
    const Eigen::VectorXd &C, double off,
    Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff, Eigen::VectorXi &Ooff) {
  // Run marching tets.
  Eigen::VectorXi I;
  marching_tets(V, TT, C, off, Voff, Foff, I);

  // Need to update the marker values.
  Ooff.resize(Voff.rows());
  for (int i = 0; i < Ooff.rows(); ++i) {
    if (I(i) == -1) {
      Ooff(i) = GLOBAL::nonoriginal_marker;
    } else {
      Ooff(i) = TO(I(i));
    }
  }

  // Dirty, dirty hack.
  Eigen::VectorXi O_new = Ooff;
  quickHack(Voff, Foff, O_new);
  //Helpers::viewTriMesh(Voff, Foff, O_new);
  //Helpers::viewTriMesh(Voff, Foff, Ooff);
  /*
  Helpers::extractManifoldPatch(Voff, Foff, Ooff, 5);
  Helpers::writeMeshWithMarkers("mesh_post_mtets", Voff, Foff, Ooff);

  Eigen::MatrixXd V_old = Voff;
  Eigen::MatrixXi F_old = Foff;
  Eigen::VectorXi O_old = Ooff;
  Eigen::MatrixXd V_new;
  Eigen::MatrixXi F_new;
  Eigen::VectorXi O_new;

  int nrem;
  int numcalls = 0;
  // First, remove all inner-ring vertices with valence 3:
  while ((nrem = removeRing(V_old, F_old, O_old, V_new, F_new, O_new, 1)) > 0) {
    numcalls++;
    printf("type 1 Removed %d\n", nrem);
    Helpers::collapseSmallTriangles(V_new, F_new);
    Helpers::removeUnreferenced(V_new, F_new, O_new);
    Helpers::viewTriMesh(V_new, F_new, O_new);
    V_old = V_new;
    F_old = F_new;
    O_old = O_new;
  }
  // Then, remove all outer-ring vertices with valence 3:
  while ((nrem = removeRing(V_old, F_old, O_old, V_new, F_new, O_new, 2)) > 0) {
    numcalls++;
    printf("type 2 Removed %d\n", nrem);
    Helpers::collapseSmallTriangles(V_new, F_new);
    Helpers::removeUnreferenced(V_new, F_new, O_new);
    Helpers::viewTriMesh(V_new, F_new, O_new);
    V_old = V_new;
    F_old = F_new;
    O_old = O_new;
  }
  // Then, remove everything else.
  while ((nrem = removeRing(V_old, F_old, O_old, V_new, F_new, O_new, 0)) > 0) {
    numcalls++;
    printf("type 0 Removed %d\n", nrem);
    Helpers::collapseSmallTriangles(V_new, F_new);
    Helpers::removeUnreferenced(V_new, F_new, O_new);
    Helpers::viewTriMesh(V_new, F_new, O_new);
    V_old = V_new;
    F_old = F_new;
    O_old = O_new;
  }

  printf("Number of times removeRing called: %d\n", numcalls);
  Voff = V_new;
  Foff = F_new;
  Ooff = O_new;
  */
  Helpers::removeUnreferenced(Voff, Foff, Ooff);
  //Helpers::viewTriMesh(Voff, Foff, Ooff);

  Eigen::VectorXi tmp;
  Helpers::extractManifoldPatch(Voff, Foff, Ooff, 5, true);
  Helpers::collapseSmallTriangles(Voff, Foff);
  Helpers::removeUnreferenced(Voff, Foff, Ooff);
  Helpers::writeMeshWithMarkers("mesh_post_rring", Voff, Foff, Ooff);

  if (!Helpers::isMeshOkay(Voff, Foff)) {
    cout << __LINE__ << endl;
    Helpers::viewTriMesh(Voff, Foff, Ooff);
    cout << "Marching offset surface mesh is ill-posed." << endl;
    exit(1);
  }
}

} // namespace OffsetSurface
