#include "OpenMeshConversion.h"

namespace Remeshing {

void iglToOpenMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
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

} // namspace Remeshing
