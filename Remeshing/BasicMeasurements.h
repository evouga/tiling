#ifndef BASIC_MEASUREMENTS_H
#define BASIC_MEASUREMENTS_H

#include "OpenMeshGlobals.h"

namespace geometry{

double squaredDistancePointToTriangle(			OpenMesh::Vec3f & _point,
												OpenMesh::Vec3f & _triangleVertexA,
												OpenMesh::Vec3f & _triangleVertexB,
												OpenMesh::Vec3f & _triangleVertexC);
												
double squaredDistancePointToTriangle(			OpenMesh::Vec3f & _point,
												OpenMesh::Vec3f & _projectedPoint,
												OpenMesh::Vec3f & _triangleVertexA,
												OpenMesh::Vec3f & _triangleVertexB,
												OpenMesh::Vec3f & _triangleVertexC);
												
double squaredDistancePointToTriangle(			OpenMesh::Vec3f & _point,
												OpenMesh::Vec3f & _projectedPoint,
												OpenMesh::Vec3f & _projectedPointBaryCentricCoordinates,
												OpenMesh::Vec3f & _triangleVertexA,
												OpenMesh::Vec3f & _triangleVertexB,
												OpenMesh::Vec3f & _triangleVertexC);
												
double signedNormalAngleOfAdjacentTriangles(	TriangleMesh::HalfedgeHandle & _he,
												TriangleMesh::FaceHandle & _f1,
												TriangleMesh::FaceHandle & _f2,
												TriangleMesh * inputMesh);

double triangleSurfaceArea(						TriangleMesh::FaceHandle & _f1,
												TriangleMesh * inputMesh);

double triangleSurfaceArea(						OpenMesh::Vec3f & pointA,
												OpenMesh::Vec3f & pointB,
												OpenMesh::Vec3f & pointC);
							

double triangleAngleBAC(						OpenMesh::Vec3f & _pointB, 
												OpenMesh::Vec3f & _pointA,
												OpenMesh::Vec3f & _pointC);

}

#endif
