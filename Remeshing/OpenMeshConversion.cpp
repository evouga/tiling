#include "OpenMeshConversion.h"

namespace Remeshing {

void IGLToOpenMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                   TriangleMesh &tMesh) {
  // Add vertices.
  std::vector<TriangleMesh::VertexHandle> vhandle;
  for (int i = 0; i < V.rows(); ++i) {
    vhandle.push_back(tMesh.add_vertex(TriangleMesh::Point(V(i, 0), V(i, 1), V(i, 2))));
    tMesh.data(vhandle[i]).setOriginalIndex(i);
  }

  // Add faces.
  for (int i = 0; i < F.rows(); ++i) {
    std::vector<TriangleMesh::VertexHandle> face_vhandles;
    for (int j = 0; j < F.cols(); ++j) {
      face_vhandles.push_back(vhandle[F(i, j)]);
    }
    tMesh.add_face(face_vhandles);
  }
}

void OpenMeshToIGL(const TriangleMesh &tMesh,
                   Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &origIdx) {
  V.resize(tMesh.n_vertices(), 3);
  F.resize(tMesh.n_faces(), 3);
  origIdx.resize(tMesh.n_vertices());

  // Look at each of the vertices.
  for (int i = 0; i < tMesh.n_vertices(); ++i) {
    TriangleMesh::VertexHandle vv = TriangleMesh::VertexHandle(i);
    const TriangleMesh::Point &pt = tMesh.point(vv);

    V.row(i) << pt[0], pt[1], pt[2];
    int orig_idx = tMesh.data(vv).getOriginalIndex();
    origIdx(i) = orig_idx;
  }

  // Navigate the mesh to find faces.
  int idx = 0;
  for (TriangleMesh::FaceIter f_it = tMesh.faces_begin();
       f_it != tMesh.faces_end();
       ++f_it, ++idx) {
    if (tMesh.valence(f_it.handle()) != 3) {
      printf("Warning: face with valence not 3 (%d)\n", 
             tMesh.valence(f_it.handle()));
    }
    TriangleMesh::ConstFaceVertexIter fvIt = tMesh.cfv_iter(*f_it);
    for (int j = 0; j < 3; ++j, ++fvIt) {
      F(idx, j) = fvIt->idx();
    }
  }
}

} // namspace Remeshing
