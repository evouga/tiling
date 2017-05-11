#include "offsetSurface.h"

#include <vector>
#include <map>
#include <utility>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <igl/boundary_facets.h>
#include <igl/copyleft/cgal/complex_to_mesh.h>
#include <igl/exterior_edges.h>
#include <igl/in_element.h>
#include <igl/unique.h>
#include <igl/collapse_edge.h>
#include <igl/unique_edge_map.h>
#include <igl/edge_flaps.h>
#include <igl/adjacency_list.h>

#include "glob_defs.h"
#include "Helpers.h"
#include "marching_tets.h"

#define OFFSET_DEBUG 0

using namespace std;

// begin namespace for helper methods.
namespace {

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
int removeRing(const Eigen::MatrixXd &V_old, const Eigen::MatrixXi &F_old,
               const Eigen::VectorXi &O_old,
               Eigen::MatrixXd &V_new, Eigen::MatrixXi &F_new,
               Eigen::VectorXi &O_new) {
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
    if (removeEdge(a, b, F_with_dead_faces, vertex_to_faces))
      have_removed.insert(b);
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

  Helpers::extractManifoldPatch(Voff, Foff, Ooff, 5);

  Eigen::MatrixXd V_old = Voff;
  Eigen::MatrixXi F_old = Foff;
  Eigen::VectorXi O_old = Ooff;
  Eigen::MatrixXd V_new;
  Eigen::MatrixXi F_new;
  Eigen::VectorXi O_new;

  while (removeRing(V_old, F_old, O_old, V_new, F_new, O_new) > 0) {
    V_old = V_new;
    F_old = F_new;
    O_old = O_new;
  }

  ///////////////////////////////////////////////////////////////////
  // Debug.
  cout << (Voff.rows() - V_new.rows()) << " vertices removed." << endl;
  cout << (Foff.rows() - F_new.rows()) << " faces removed." << endl;
  ///////////////////////////////////////////////////////////////////

  Voff = V_new;
  Foff = F_new;
  Ooff = O_new;

  ///////////////////////////////////////////////////////////////////
  // Debug.
  Eigen::VectorXi tmp;
  cout << "manifold before:" << Helpers::is_edge_manifold(Foff, tmp) << endl;
  Helpers::extractManifoldPatch(Voff, Foff, Ooff, 5, false);
  cout << "manifold after:" << Helpers::is_edge_manifold(Foff, tmp) << endl;
  Helpers::extractManifoldPatch(Voff, Foff, Ooff, 5, false);
  ///////////////////////////////////////////////////////////////////

  Helpers::collapseSmallTriangles(Voff, Foff);

  cout << "debugging" << endl;
  cout << "is okay: " << Helpers::isMeshOkay(Voff, Foff) << endl;

  Helpers::removeUnreferenced(Voff, Foff, Ooff);

  cout << "is okay: " << Helpers::isMeshOkay(Voff, Foff) << endl;

  // Helpers::viewTriMesh(Voff, Foff, Ooff);
}

} // namespace OffsetSurface
