#include "Remesh.h"

#include "../glob_defs.h"
#include "IsotropicRemeshing.h"
#include "OpenMeshConversion.h"

namespace {
double avg_edgeLength(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const Eigen::VectorXi &M) {
  double total_el = 0;
  int num_counted = 0;
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < F.cols(); ++j) {
      if (M(F(i, j)) == GLOBAL::nonoriginal_marker) continue;

      num_counted++;
      int v1 = F(i, j);
      int v2 = F(i, (j+1) % F.cols());
      total_el += (V.row(v1) - V.row(v2)).norm();
    }
  }

  return total_el / num_counted;
}
} // namespace

namespace Remeshing {
void remesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &M,
            Eigen::MatrixXd &Vout, Eigen::MatrixXi &Fout, Eigen::VectorXi &Mout,
            double e_length) {
  // The boundary vertices.
  double zmin = 1e10;
  double zmax = -1e10;
  for (int i = 0; i < V.rows(); i++) {
    zmin = std::min(zmin, V(i, 2));
    zmax = std::max(zmax, V(i, 2));
  }

  // Compute average edge length.
  if (e_length <= 0) {
    e_length = avg_edgeLength(V, F, M);
  }

  // Convert to OpenMesh.
  TriangleMesh tri;
  IGLToOpenMesh(V, F, tri);
  // Do some more setup.
  tri.update_face_normals();
  tri.update_vertex_normals();
  tri.request_vertex_status();
  tri.request_edge_status();
  tri.request_face_status();

  // Check all edges to see if they're protected or not.
  TriangleMesh::EIter e_it;
  for (e_it = tri.edges_begin(); e_it != tri.edges_end(); ++e_it) {
    // Get both vertices.
    TriangleMesh::VHandle v0 = tri.to_vertex_handle(tri.halfedge_handle(*e_it, 0));
    TriangleMesh::VHandle v1 = tri.to_vertex_handle(tri.halfedge_handle(*e_it, 1));
    // See if both are nonoriginal
    if (M(tri.data(v0).getOriginalIndex()) != GLOBAL::nonoriginal_marker &&
        M(tri.data(v1).getOriginalIndex()) != GLOBAL::nonoriginal_marker) {
      tri.data(*e_it).setProtected(true);
    }
    // Also set vertices if they're original.
    if (M(tri.data(v0).getOriginalIndex()) != GLOBAL::nonoriginal_marker) {
      tri.data(v0).setProtected(true);
    }
    if (M(tri.data(v1).getOriginalIndex()) != GLOBAL::nonoriginal_marker) {
      tri.data(v1).setProtected(true);
    }
  }
  
  /*
  // Check for bad vertices.
  for (int i = 0; i < tri.n_vertices(); ++i) {
    TriangleMesh::VertexHandle vv = TriangleMesh::VertexHandle(i);
    const TriangleMesh::Point &pt = tMesh.point(vv);

    // Set statitonary if original marker set.
    if (M(tri.data(vv).getOriginalIndex()) == GLOBAL::original_marker) {
      
    }
  }
  */

  // Then, do the remeshing--but only split long edge, collapse short edges, 
  // and equalize valences.
  IsotropicRemeshing ir(e_length);
  int types = ISOTROPIC_REMESHING_TYPES::SPLIT_LONG_EDGES |
              ISOTROPIC_REMESHING_TYPES::COLLAPSE_SHORT_EDGES |
              ISOTROPIC_REMESHING_TYPES::EQUALIZE_VALENCES;
  ir.remesh(&tri, 1, types);

  // Then, copy it all back.
  Eigen::VectorXi origIdx;
  OpenMeshToIGL(tri, Vout, Fout, origIdx);
  // Set the new markers as well.
  Mout.resize(Vout.rows());
  for (int i = 0; i < Vout.rows(); ++i) {
    // If it doesn't correspond to a previous vertex, it's a "nonoriginal".
    if (origIdx(i) == -1) {
      Mout(i) = GLOBAL::nonoriginal_marker;
    } else {
      // Otherwise, it's the same marker as it was before.
      Mout(i) = M(origIdx(i));
      int oidx = origIdx(i);
      double moved = ( V.row(oidx) - Vout.row(i) ).norm();
      if ( moved > GLOBAL::EPS && M(oidx) != GLOBAL::nonoriginal_marker) {
        printf("point %d(m:%d new:%d) moved by %lf\n", oidx, M(oidx), i, moved);
        Mout(i) = 10;
      }
      if (V(i, 2) < zmin || V(i, 2) > zmax) {
        printf("point %d(%d) moved out of bounds! by %lf vs %lf,%lf\n",
               oidx, M(oidx), V(i, 2), zmin, zmax);
      }
    }
  }
}

} // namespace Remeshing
