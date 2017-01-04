#include "offsetSurface.h"

#include <vector>
#include <map>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <igl/boundary_facets.h>
#include <igl/copyleft/cgal/complex_to_mesh.h>
#include <igl/exterior_edges.h>
#include <igl/in_element.h>

#include "glob_defs.h"

#define OFFSET_DEBUG 0

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

void generateOffsetSurface_naive(const Eigen::MatrixXd &V,
                                 const Eigen::MatrixXi &TT,
                                 const Eigen::VectorXi &TO,
                                 const Eigen::VectorXd &C, double off,
                                 Eigen::MatrixXd &Voff, Eigen::MatrixXi &Foff, Eigen::VectorXi &Ooff) {
  std::vector<int> s;
  for (int i = 0; i < TT.rows(); ++i) {
    if (C(TT(i, 0)) <= off && C(TT(i, 1)) <= off &&
        C(TT(i, 2)) <= off && C(TT(i, 3)) <= off) {
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
  // And remove duplicate faces.
  Eigen::MatrixXi Fu;
  igl::unique_rows(Foff, Fu,unique_to_all,all_to_unique);
  Foff = Fu;
}

} // namespace OffsetSurface
