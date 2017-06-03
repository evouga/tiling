#ifndef NEAREST_TRIANGLE_SEARCH_H
#define NEAREST_TRIANGLE_SEARCH_H

#include "OpenMeshGlobals.h"

#include "KDTreeNode.h"
#include "BasicMeasurements.h"

#include <vector>

namespace dataStructure{
class NearestTriangleSearch{

	public:
		NearestTriangleSearch();

		void initializeKDTreeBasedSearchStructure(TriangleMesh * inputMesh,unsigned int maxDepth, unsigned int minNumberOfTriangles);
		void destroyKDTreeBasedSearchStructure();

		void computeClosestTriangleOfPoint(OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle);
		
	private:	

		void destroyNode(KDTreeNode * node);

		//void getFarthestAncestorNodeWithPossibleNearestTriangleOfPoint(KDTreeNode * nodeOfSubSpace,KDTreeNode ** farthestAncestorNode,OpenMesh::Vec3f & queryPoint, double minDistanceOfAllTriangleInLeafNode,unsigned int counter);

		void closestTriangleOfPointInAncestorCells(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,double * minDistance,TriangleMesh::FaceHandle & closestTriangle);
		
		void closestTriangleOfPointInDescendantCells(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle, double minDistanceOfAllTriangleInLeafNode, double * minDistance);

		double closestTriangleOfPointInCell(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle);

		void buildKDTreeWithMedianSplit(TriangleMesh * inputMesh,unsigned int maxDepth, unsigned int minNumberOfTriangles);

		void subdivideSpace(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace,unsigned int maxDepth, unsigned int minNumberOfTriangles,unsigned int currentKDTreeDepth);

		void computeCoordinateOfSplitting(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace);

		bool terminateKDTreeConstruction(KDTreeNode * nodeOfSubSpace, unsigned int currentKDTreeDepth, unsigned int maxDepth, unsigned int minNumberOfTriangles);

		void computeCellBoundingBox(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & minPoint, OpenMesh::Vec3f & maxPoint); // currently working on

		void sortTrianglesIntoChildCells(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace);

		KDTreeNode * findLeafWithQueryPoint(OpenMesh::Vec3f & queryPoint,KDTreeNode * entryNode);


		KDTreeNode * kDTreeRootNode;
	

};
}

#endif
