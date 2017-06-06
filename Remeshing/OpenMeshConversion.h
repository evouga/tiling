#ifndef REMESHING_OPEN_MESH_CONVERSION_H
#define REMESHING_OPEN_MESH_CONVERSION_H

#include <Eigen/Core>

#include "MeshTraits.h"

namespace Remeshing {
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriangleMesh;                

void IGLToOpenMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                   TriangleMesh &tMesh);
void OpenMeshToIGL(const TriangleMesh &tMesh,
                   Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &origIdx);

} // namespace OpenMesh
#endif // REMESHING_OPEN_MESH_CONVERSION_H
