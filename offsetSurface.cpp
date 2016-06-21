#include "offsetSurface.h"

#include <vector>

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>

#include <igl/copyleft/cgal/complex_to_mesh.h>
#include <igl/in_element.h>

#define HIGHEST_TEMP  100.

#define OFFSET_DEBUG 1
//#define OFF_MINUS

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

  // Code borrowed from igl: https://github.com/libigl/libigl/blob/master/include/igl/copyleft/cgal/signed_distance_isosurface.cpp
  const auto &Vmax = V.colwise().maxCoeff();
  const auto &Vmin = V.colwise().minCoeff();
  const double bbd = (Vmax-Vmin).norm();
  const double r = bbd/2.;
  const auto &Vmid = 0.5*(Vmax + Vmin);

  Point_3 cmid(Vmid(0), Vmid(1), Vmid(2));
  // Need the triangulation and the offset
  Function fun = 
      [&](const Point_3 &qq) -> FT
      {
        Eigen::RowVectorXd q(3);
        q << qq.x(), qq.y(), qq.z();
        // Only return first element containing q
        std::vector<int> tets = tree.find(V, TT, q, true);
        if (tets.size() < 1) {
          return HIGHEST_TEMP;
        }
        int tet = tets[0];
        // Make sure it's in a tet

        // Get the 4 vertices
        Eigen::MatrixXd mat(4, 4);
        Eigen::Vector4d query(4);
        for (int i = 0; i < 4; ++i) {
          mat(0, i) = V(TT(tet,i), 0);
          mat(1, i) = V(TT(tet,i), 1);
          mat(2, i) = V(TT(tet,i), 2);
          mat(3, i) = 1;
          if (OFFSET_DEBUG) 
            printf("%f %f %f %f 1\n",
                   V(TT(tet,i), 0), V(TT(tet,i), 1), V(TT(tet,i), 2),
                   C(TT(tet,i)));
        }
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
             << alpha << "\n";
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
#ifdef OFF_MINUS
        return off - val;
#else
        return val - off;
#endif
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
  // Now, find the surface
  //Sphere_3 bounding_sphere(cmid, r*r);  // squared radius
  Sphere_3 bounding_sphere(cmid, (r+off)*(r+off));  // squared radius
  Surface_3 surface(fun, bounding_sphere);
  // Use default values for angle bound, radius bound, distance bound (respectively)
  CGAL::Surface_mesh_default_criteria_3<CTr> criteria(28, 0.2, 0.01);
  //CGAL::Surface_mesh_default_criteria_3<CTr> criteria(28, 1.1, 1.1);
  // Mesh surface
  CTr tr;  // 3D-Delaunay triangulation
  C2t3 c2t3(tr); // 2D-complex in 3D-Delaunay triangulation
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  // Complex to (V,F)
  igl::copyleft::cgal::complex_to_mesh(c2t3, Voff, Foff);
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

}
