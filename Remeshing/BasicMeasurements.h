#ifndef BASIC_MEASUREMENTS_H
#define BASIC_MEASUREMENTS_H

#include "OpenMeshGlobals.h"

namespace geometry{

double squaredDistancePointToTriangle(			TriangleMesh::Point & _point,
												TriangleMesh::Point & _triangleVertexA,
												TriangleMesh::Point & _triangleVertexB,
												TriangleMesh::Point & _triangleVertexC);
												
double squaredDistancePointToTriangle(			TriangleMesh::Point & _point,
												TriangleMesh::Point & _projectedPoint,
												TriangleMesh::Point & _triangleVertexA,
												TriangleMesh::Point & _triangleVertexB,
												TriangleMesh::Point & _triangleVertexC);
												
double squaredDistancePointToTriangle(			TriangleMesh::Point & _point,
												TriangleMesh::Point & _projectedPoint,
												TriangleMesh::Point & _projectedPointBaryCentricCoordinates,
												TriangleMesh::Point & _triangleVertexA,
												TriangleMesh::Point & _triangleVertexB,
												TriangleMesh::Point & _triangleVertexC);
												
double signedNormalAngleOfAdjacentTriangles(	TriangleMesh::HalfedgeHandle & _he,
												TriangleMesh::FaceHandle & _f1,
												TriangleMesh::FaceHandle & _f2,
												TriangleMesh * inputMesh);

double triangleSurfaceArea(						TriangleMesh::FaceHandle & _f1,
												TriangleMesh * inputMesh);

double triangleSurfaceArea(						TriangleMesh::Point & pointA,
												TriangleMesh::Point & pointB,
												TriangleMesh::Point & pointC);
							

double triangleAngleBAC(						TriangleMesh::Point & _pointB, 
												TriangleMesh::Point & _pointA,
												TriangleMesh::Point & _pointC);

}

#endif
