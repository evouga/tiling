#ifndef OPEN_MESH_GLOBALS_H
#define OPEN_MESH_GLOBALS_H

#undef check

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>

#include "MeshTraits.h"

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriangleMesh;

#endif //OPEN_MESH_GLOBALS_H 
