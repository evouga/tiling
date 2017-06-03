#ifndef REMESHING_OPEN_MESH_CONVERSION_H
#define REMESHING_OPEN_MESH_CONVERSION_H

#include <Eigen/Core>

#include "MeshTraits.h"

namespace Remeshing {
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriangleMesh;                

void iglToOpenMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXd &F,
                   TriangleMesh &tMesh);

} // namespace OpenMesh
#endif // REMESHING_OPEN_MESH_CONVERSION_H
