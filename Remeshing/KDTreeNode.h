#ifndef KD_TREE_NODE_H
#define KD_TREE_NODE_H

#include "OpenMeshGlobals.h"

#include <vector>

namespace dataStructure{

typedef struct KDTreeAxiAlignedBoundingBox_{
		
  TriangleMesh::Point minCorner;
  TriangleMesh::Point maxCorner;

}KDTreeAxiAlignedBoundingBox;


enum COORDINATE_AXIS{
	X_AXIS=0,
	Y_AXIS,
	Z_AXIS,
};

typedef struct KDTreeNode_{
		
	bool isLeaf;

	KDTreeAxiAlignedBoundingBox boundingBox;

	enum COORDINATE_AXIS coordinateAxisOfSplitting;

	double splittingPlaneCoordinate;

	std::vector<TriangleMesh::FaceHandle> * faceHandles;

	struct KDTreeNode_ * leftChild;
	struct KDTreeNode_ * rightChild;
	struct KDTreeNode_ * parent;
	
} KDTreeNode;
}
#endif
