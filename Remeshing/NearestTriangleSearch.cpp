#include "NearestTriangleSearch.h"
#include "BasicMeasurements.h"

#include <algorithm> 
#include <limits>
namespace dataStructure{

NearestTriangleSearch::NearestTriangleSearch(){


}

void NearestTriangleSearch::initializeKDTreeBasedSearchStructure(TriangleMesh * inputMesh,unsigned int maxDepth,  unsigned int minNumberOfTriangles){
		
	buildKDTreeWithMedianSplit(inputMesh,maxDepth,minNumberOfTriangles);
	
}

void NearestTriangleSearch::destroyKDTreeBasedSearchStructure(){
	destroyNode(kDTreeRootNode);
	delete kDTreeRootNode;
	kDTreeRootNode = NULL;
}

void NearestTriangleSearch::destroyNode(KDTreeNode * node){
	if(node->isLeaf){
		if (node->faceHandles!=NULL){
			node->faceHandles->clear();
			delete node->faceHandles;
			node->faceHandles=NULL;
		}
	}
	else if(node->leftChild->isLeaf && node->rightChild->isLeaf){
		destroyNode(node->leftChild);
		destroyNode(node->rightChild);
		delete node->leftChild;
		delete node->rightChild;
		node->leftChild=NULL;
		node->rightChild=NULL;
		node->isLeaf = true;
		destroyNode(node);
	}
	else{
		destroyNode(node->leftChild);
		destroyNode(node->rightChild);
	}
}


void NearestTriangleSearch::computeClosestTriangleOfPoint(OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle){
	/*
	KDTreeNode * leafWithQueryPoint;
	KDTreeNode * farthestAncestorNode = new KDTreeNode;
	
	// recursively traverse kdtree to find leaf containing the query point
	
	leafWithQueryPoint = findLeafWithQueryPoint(queryPoint,kDTreeRootNode);

	// find minimum distance of all triangles in leaf node
	
	double minDistanceOfAllTriangleInLeafNode = closestTriangleOfPointInCell(leafWithQueryPoint,queryPoint,inputMesh,closestTriangle);
	
	//track back up the tree to see if an ancestor cell covers another possible nearest triangle
	farthestAncestorNode = leafWithQueryPoint;
	unsigned int counter = 0;
	//printf("=======\n");
	getFarthestAncestorNodeWithPossibleNearestTriangleOfPoint(leafWithQueryPoint,&farthestAncestorNode,queryPoint,minDistanceOfAllTriangleInLeafNode,counter);
	double temp = std::numeric_limits<float>::max();
	closestTriangleOfPointInDescendantCells(farthestAncestorNode,queryPoint,inputMesh,closestTriangle,minDistanceOfAllTriangleInLeafNode,&temp);
	 */
	
	KDTreeNode * leafWithQueryPoint = findLeafWithQueryPoint(queryPoint,kDTreeRootNode);
	double minDistanceOfAllTriangleInLeafNode = closestTriangleOfPointInCell(leafWithQueryPoint,queryPoint,inputMesh,closestTriangle);
	closestTriangleOfPointInAncestorCells(leafWithQueryPoint,queryPoint,inputMesh,&minDistanceOfAllTriangleInLeafNode,closestTriangle);
	
}
/*
void NearestTriangleSearch::getFarthestAncestorNodeWithPossibleNearestTriangleOfPoint(KDTreeNode * nodeOfSubSpace,KDTreeNode ** farthestAncestorNode,OpenMesh::Vec3f & queryPoint, double minDistanceOfAllTriangleInLeafNode, unsigned int counter){
	
	if(nodeOfSubSpace->parent != NULL){
		if(fabs(nodeOfSubSpace->parent->splittingPlaneCoordinate-queryPoint[nodeOfSubSpace->parent->coordinateAxisOfSplitting])<minDistanceOfAllTriangleInLeafNode){
			
			*farthestAncestorNode = nodeOfSubSpace->parent;
			
			//if (counter>5){
			//printf("counter: %d\n",counter);
			//printf("mindistanceToTriangleInLeafNode: %f, ",minDistanceOfAllTriangleInLeafNode);
			//printf("splittingPlaneCoordinate: %f, ",nodeOfSubSpace->parent->splittingPlaneCoordinate);
			//printf("queryPoint: %f\n",queryPoint[nodeOfSubSpace->parent->coordinateAxisOfSplitting]);
			//}
		}
		getFarthestAncestorNodeWithPossibleNearestTriangleOfPoint(nodeOfSubSpace->parent,farthestAncestorNode,queryPoint,minDistanceOfAllTriangleInLeafNode,counter+1);
	}
	
}
*/

void NearestTriangleSearch::closestTriangleOfPointInAncestorCells(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,double * minDistance,TriangleMesh::FaceHandle & closestTriangle)
{
	if(nodeOfSubSpace->parent != NULL){
		if(fabs(nodeOfSubSpace->parent->splittingPlaneCoordinate-queryPoint[nodeOfSubSpace->parent->coordinateAxisOfSplitting])<*minDistance){
			double newMinDistance = std::numeric_limits<float>::max();
			TriangleMesh::FaceHandle newClosestTriangle;
			if(nodeOfSubSpace==nodeOfSubSpace->parent->leftChild){
				closestTriangleOfPointInDescendantCells(nodeOfSubSpace->parent->rightChild,queryPoint,inputMesh,newClosestTriangle,*minDistance,&newMinDistance);
			}
			else
			{
				closestTriangleOfPointInDescendantCells(nodeOfSubSpace->parent->leftChild,queryPoint,inputMesh,newClosestTriangle,*minDistance,&newMinDistance);
			}
			if(newMinDistance<*minDistance)
			{
				*minDistance = newMinDistance;
				closestTriangle = newClosestTriangle;
			}
			
		}
		closestTriangleOfPointInAncestorCells(nodeOfSubSpace->parent,queryPoint,inputMesh,minDistance,closestTriangle);
	}

}

void NearestTriangleSearch::closestTriangleOfPointInDescendantCells(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle, double minDistanceOfAllTriangleInLeafNode, double * minDistance){

	if(!nodeOfSubSpace->isLeaf){
		double minDistanceInLeftCell =  std::numeric_limits<float>::max();
		double minDistanceInRightCell =  std::numeric_limits<float>::max();;
		
		TriangleMesh::FaceHandle closestTriangleInLeftCell;
		TriangleMesh::FaceHandle closestTriangleInRightCell;
		
		// if splitting plane intersects ball of query point with radius the smallest distance to leaf
		if(fabs(nodeOfSubSpace->splittingPlaneCoordinate-queryPoint[nodeOfSubSpace->coordinateAxisOfSplitting])<minDistanceOfAllTriangleInLeafNode){
			
			closestTriangleOfPointInDescendantCells(nodeOfSubSpace->leftChild,queryPoint,inputMesh,closestTriangleInLeftCell, minDistanceOfAllTriangleInLeafNode,&minDistanceInLeftCell);
			closestTriangleOfPointInDescendantCells(nodeOfSubSpace->rightChild,queryPoint,inputMesh,closestTriangleInRightCell, minDistanceOfAllTriangleInLeafNode,&minDistanceInRightCell);
			
			if(minDistanceInLeftCell<minDistanceInRightCell){
				*minDistance = minDistanceInLeftCell;
				closestTriangle = closestTriangleInLeftCell;
			}
			else{
				*minDistance = minDistanceInRightCell;
				closestTriangle = closestTriangleInRightCell;
			}
			
		}
		else if(queryPoint[nodeOfSubSpace->coordinateAxisOfSplitting] < nodeOfSubSpace->splittingPlaneCoordinate)
		{
			closestTriangleOfPointInDescendantCells(nodeOfSubSpace->leftChild,queryPoint,inputMesh,closestTriangle, minDistanceOfAllTriangleInLeafNode,minDistance);
		}
		else if(queryPoint[nodeOfSubSpace->coordinateAxisOfSplitting] > nodeOfSubSpace->splittingPlaneCoordinate)
		{
			closestTriangleOfPointInDescendantCells(nodeOfSubSpace->rightChild,queryPoint,inputMesh,closestTriangle, minDistanceOfAllTriangleInLeafNode,minDistance);
		}
	}
	else{
		*minDistance = closestTriangleOfPointInCell(nodeOfSubSpace,queryPoint,inputMesh,closestTriangle);
	}
}

double NearestTriangleSearch::closestTriangleOfPointInCell(KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & queryPoint,TriangleMesh * inputMesh,TriangleMesh::FaceHandle & closestTriangle){
	TriangleMesh::ConstFaceVertexIter cfvIt;
	double minQueryPointTriangleDistance = std::numeric_limits<double>::max();
	TriangleMesh::FaceHandle nearestTriangleHandle;
	
	for(unsigned int i=0;i<nodeOfSubSpace->faceHandles->size();i++){
		
		TriangleMesh::FaceHandle currentTriangleHandle = nodeOfSubSpace->faceHandles->at(i);

		cfvIt = inputMesh->cfv_iter(currentTriangleHandle);

		TriangleMesh::Point vertex1 = inputMesh->point( *cfvIt );
		++cfvIt;
		TriangleMesh::Point vertex2 = inputMesh->point( *cfvIt );
		++cfvIt;
		TriangleMesh::Point vertex3 = inputMesh->point( *cfvIt );

		double currentDistance = sqrt(geometry::squaredDistancePointToTriangle(queryPoint,vertex1,vertex2,vertex3));
		
		if(currentDistance<minQueryPointTriangleDistance){
				nearestTriangleHandle = currentTriangleHandle;
				minQueryPointTriangleDistance = currentDistance;
		}		
	}
	closestTriangle = nearestTriangleHandle;
	return minQueryPointTriangleDistance;
} 

KDTreeNode * NearestTriangleSearch::findLeafWithQueryPoint(OpenMesh::Vec3f & queryPoint,KDTreeNode * entryNode){
	if(!entryNode->isLeaf){
		
		enum COORDINATE_AXIS currentCoordinateAxisOfSplitting = entryNode->coordinateAxisOfSplitting;
		double currentSplittingPlaneCoordinate = entryNode->splittingPlaneCoordinate;

		if(queryPoint[currentCoordinateAxisOfSplitting]<currentSplittingPlaneCoordinate){
			return findLeafWithQueryPoint(queryPoint,entryNode->leftChild);
		}
		else{
			return findLeafWithQueryPoint(queryPoint,entryNode->rightChild);
		}
	}
	else{
		//printf("entry node number of triangles %d\n",entryNode->faceHandles->size());
		return entryNode;
	}
}

void NearestTriangleSearch::buildKDTreeWithMedianSplit(TriangleMesh * inputMesh,unsigned int maxDepth, unsigned int minNumberOfTriangles){
	
		unsigned int currentKDTreeDepth = 0;

		// initialize root node
		kDTreeRootNode = new KDTreeNode;
		kDTreeRootNode->isLeaf = true;
		kDTreeRootNode->parent = NULL;
		kDTreeRootNode->coordinateAxisOfSplitting = X_AXIS;
		kDTreeRootNode->splittingPlaneCoordinate = 0.0;
		
		TriangleMesh::FaceIter fIt=inputMesh->faces_begin();
		TriangleMesh::FaceIter fEnd=inputMesh->faces_end();
		kDTreeRootNode->faceHandles = new std::vector<TriangleMesh::FaceHandle> ;
		kDTreeRootNode->faceHandles->clear();
		for(;fIt!=fEnd;++fIt){
			kDTreeRootNode->faceHandles->push_back( *fIt );
		}

		subdivideSpace(inputMesh,kDTreeRootNode,maxDepth,minNumberOfTriangles,currentKDTreeDepth);
}

void NearestTriangleSearch::subdivideSpace(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace,unsigned int maxDepth, unsigned int minNumberOfTriangles,unsigned int currentKDTreeDepth){

	if(!terminateKDTreeConstruction(nodeOfSubSpace,currentKDTreeDepth,maxDepth,minNumberOfTriangles)){
	//printf("current depth: %d, number of triangles: %d",currentKDTreeDepth,nodeOfSubSpace->faceHandles->size());
			
		// subdivide node //////////////////////////////////////////////////////
		
		// update global state
		// currentKDTreeDepth++;
		
		// update current Node	
		nodeOfSubSpace->isLeaf = false;
		
		// compute coordinateAxis of splitting
		OpenMesh::Vec3f cellBoundingBoxMin;
		OpenMesh::Vec3f cellBoundingBoxMax;
		
		computeCellBoundingBox(inputMesh,nodeOfSubSpace,cellBoundingBoxMin,cellBoundingBoxMax);
		double lengthX, lengthY, lengthZ;
		
			// find longest bounding box direction in xyz-axes
			lengthX = fabs(cellBoundingBoxMax[0]-cellBoundingBoxMin[0]);
			lengthY = fabs(cellBoundingBoxMax[1]-cellBoundingBoxMin[1]);
			lengthZ = fabs(cellBoundingBoxMax[2]-cellBoundingBoxMin[2]);	
			//printf(", length x: %f, length y: %f, length z: %f",lengthX,lengthY,lengthZ);

			if(lengthY > lengthX)
			{
				if(lengthY > lengthZ)
				{
					nodeOfSubSpace->coordinateAxisOfSplitting=Y_AXIS;
					//printf(", Y between %f and %f",cellBoundingBoxMin[1],cellBoundingBoxMax[1]);
				}
				else
				{
					nodeOfSubSpace->coordinateAxisOfSplitting=Z_AXIS;
					//printf(", Z between %f and %f",cellBoundingBoxMin[2],cellBoundingBoxMax[2]);
				}
			}
			else
			{
				if(lengthX > lengthZ)
				{
					nodeOfSubSpace->coordinateAxisOfSplitting=X_AXIS;
					//printf(", X between %f and %f",cellBoundingBoxMin[0],cellBoundingBoxMax[0]);
				}
				else
				{
					nodeOfSubSpace->coordinateAxisOfSplitting=Z_AXIS;
					//printf(", Z between %f and %f",cellBoundingBoxMin[2],cellBoundingBoxMax[2]);
				}
			}
		
		// compute coordinate of splitting
		computeCoordinateOfSplitting(inputMesh,nodeOfSubSpace);
		//printf(", split at %f\n",nodeOfSubSpace->splittingPlaneCoordinate);

		// create children ////////////////////////////////////////////////////
		
		nodeOfSubSpace->leftChild	= new KDTreeNode;
		nodeOfSubSpace->rightChild	= new KDTreeNode;
		
		nodeOfSubSpace->leftChild->parent	= nodeOfSubSpace;
		nodeOfSubSpace->rightChild->parent	= nodeOfSubSpace;

		nodeOfSubSpace->leftChild->isLeaf	= true;
		nodeOfSubSpace->rightChild->isLeaf	= true;
		
		// cyclic next coordinate axis of splitting: slower
		/*
		if(nodeOfSubSpace->coordinateAxisOfSplitting==X_AXIS){
			nodeOfSubSpace->leftChild->coordinateAxisOfSplitting=Y_AXIS;
			nodeOfSubSpace->rightChild->coordinateAxisOfSplitting=Y_AXIS;
		
		}
		else if(nodeOfSubSpace->coordinateAxisOfSplitting==Y_AXIS){
			nodeOfSubSpace->leftChild->coordinateAxisOfSplitting=Z_AXIS;
			nodeOfSubSpace->rightChild->coordinateAxisOfSplitting=Z_AXIS;			
		}

		else if(nodeOfSubSpace->coordinateAxisOfSplitting==Z_AXIS){
			nodeOfSubSpace->leftChild->coordinateAxisOfSplitting=X_AXIS;
			nodeOfSubSpace->rightChild->coordinateAxisOfSplitting=X_AXIS;
		}
		*/
		
		// Reference Triangles into new leafs
		sortTrianglesIntoChildCells(inputMesh,nodeOfSubSpace); //TODO
		
		// free unnecessary current node memory, which are the references to the triangles
		nodeOfSubSpace->faceHandles->clear(); 
		delete nodeOfSubSpace->faceHandles;
		nodeOfSubSpace->faceHandles=NULL;
		
		// recursion /////////////////////////////////////////////////////////
			
		subdivideSpace(inputMesh,nodeOfSubSpace->leftChild,maxDepth,minNumberOfTriangles,currentKDTreeDepth+1);
		subdivideSpace(inputMesh,nodeOfSubSpace->rightChild,maxDepth,minNumberOfTriangles,currentKDTreeDepth+1);
	}
	
}


// TODO rename to cell bounding box
void NearestTriangleSearch::computeCellBoundingBox(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace,OpenMesh::Vec3f & minPoint, OpenMesh::Vec3f & maxPoint){
	TriangleMesh::FaceVertexIter fvIt;
	OpenMesh::Vec3f point0,point1,point2;
	
	OpenMesh::Vec3f currentMinPoint;
	OpenMesh::Vec3f currentMaxPoint;
	
	currentMinPoint.vectorize( std::numeric_limits< float >::max() );
	currentMaxPoint.vectorize( - std::numeric_limits< float >::max() );
	  
	for(unsigned int i=0;i<nodeOfSubSpace->faceHandles->size();i++){
		fvIt = inputMesh->fv_iter(nodeOfSubSpace->faceHandles->at(i));
		
		point0 = inputMesh->point( *fvIt );
		currentMinPoint.minimize(point0);
		currentMaxPoint.maximize(point0);

		point1 = inputMesh->point( *(++fvIt) );
		
		currentMinPoint.minimize(point1);
		currentMaxPoint.maximize(point1);
		
		point2 = inputMesh->point( *(++fvIt) );
		
		currentMinPoint.minimize(point2);
		currentMaxPoint.maximize(point2);
		
	}
	
	minPoint = currentMinPoint;
	maxPoint = currentMaxPoint;
	
}

void NearestTriangleSearch::sortTrianglesIntoChildCells(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace){
	TriangleMesh::FaceVertexIter fvIt;
	OpenMesh::Vec3f point0,point1,point2;

	double triangleLowerBound=0.0;
	double triangleUpperBound=0.0;
	
	nodeOfSubSpace->leftChild->faceHandles = new std::vector<TriangleMesh::FaceHandle>;
	nodeOfSubSpace->rightChild->faceHandles = new std::vector<TriangleMesh::FaceHandle>;
	
	nodeOfSubSpace->leftChild->faceHandles->clear();
	nodeOfSubSpace->rightChild->faceHandles->clear();

	for(unsigned int i=0;i<nodeOfSubSpace->faceHandles->size();i++){
		fvIt = inputMesh->fv_iter(nodeOfSubSpace->faceHandles->at(i));
		triangleLowerBound = std::numeric_limits<double>::max();
		triangleUpperBound = -std::numeric_limits<double>::max();
		
		point0 = inputMesh->point( *fvIt );

		if(point0[nodeOfSubSpace->coordinateAxisOfSplitting]<triangleLowerBound){
			triangleLowerBound = point0[nodeOfSubSpace->coordinateAxisOfSplitting];
		} 
		if(point0[nodeOfSubSpace->coordinateAxisOfSplitting]>triangleUpperBound){
			triangleUpperBound = point0[nodeOfSubSpace->coordinateAxisOfSplitting];
		}
		point1 = inputMesh->point( *(++fvIt) );

		if(point1[nodeOfSubSpace->coordinateAxisOfSplitting]<triangleLowerBound){
			triangleLowerBound = point1[nodeOfSubSpace->coordinateAxisOfSplitting];
		} 
		if(point1[nodeOfSubSpace->coordinateAxisOfSplitting]>triangleUpperBound){
			triangleUpperBound = point1[nodeOfSubSpace->coordinateAxisOfSplitting];
		}
		point2 = inputMesh->point( *(++fvIt) );

		if(point2[nodeOfSubSpace->coordinateAxisOfSplitting]<triangleLowerBound){
			triangleLowerBound = point2[nodeOfSubSpace->coordinateAxisOfSplitting];
		} 
		if(point2[nodeOfSubSpace->coordinateAxisOfSplitting]>triangleUpperBound){
			triangleUpperBound = point2[nodeOfSubSpace->coordinateAxisOfSplitting];
		}
		
		// sort into left children
		
		if(triangleUpperBound<nodeOfSubSpace->splittingPlaneCoordinate){
			nodeOfSubSpace->leftChild->faceHandles->push_back(nodeOfSubSpace->faceHandles->at(i));
		}	
		// sort into right children
		else if(triangleLowerBound>nodeOfSubSpace->splittingPlaneCoordinate){
			nodeOfSubSpace->rightChild->faceHandles->push_back(nodeOfSubSpace->faceHandles->at(i));	
		}
		
		// sort into left and right children
		else if((triangleLowerBound<=nodeOfSubSpace->splittingPlaneCoordinate) && (triangleUpperBound>=nodeOfSubSpace->splittingPlaneCoordinate)){
			nodeOfSubSpace->leftChild->faceHandles->push_back(nodeOfSubSpace->faceHandles->at(i));	
			nodeOfSubSpace->rightChild->faceHandles->push_back(nodeOfSubSpace->faceHandles->at(i));
		}
		
	}
	
}


void NearestTriangleSearch::computeCoordinateOfSplitting(TriangleMesh * inputMesh,KDTreeNode * nodeOfSubSpace){
	TriangleMesh::FaceVertexIter fvIt;
	
	double point0Coordinate, point1Coordinate, point2Coordinate, centerOfGravityCoordinate;
	
	std::vector<double> coordinateOfSplittingList;
	double coordinateOfSplitting=0.0;
	coordinateOfSplittingList.clear();
	enum COORDINATE_AXIS currentCoordinateAxisOfSplitting = nodeOfSubSpace->coordinateAxisOfSplitting;
	//printf("*************\n");
	for(unsigned int i=0;i<nodeOfSubSpace->faceHandles->size();i++){
		fvIt = inputMesh->fv_iter(nodeOfSubSpace->faceHandles->at(i));
		point0Coordinate = inputMesh->point( *fvIt )[currentCoordinateAxisOfSplitting];
		point1Coordinate = inputMesh->point( *(++fvIt) )[currentCoordinateAxisOfSplitting];
		point2Coordinate = inputMesh->point( *(++fvIt) )[currentCoordinateAxisOfSplitting];
		centerOfGravityCoordinate = (point0Coordinate+point1Coordinate+point2Coordinate)/3.0;
		
		coordinateOfSplittingList.push_back(centerOfGravityCoordinate);
		//printf("center of gravity: %f\n",centerOfGravityCoordinate);
	}
	sort(coordinateOfSplittingList.begin(), coordinateOfSplittingList.end()); //introsert, not quicksort :-)
	
	coordinateOfSplitting = coordinateOfSplittingList[coordinateOfSplittingList.size()/2];
	/*
	printf("*************\n");
	for(unsigned int i=0;i<coordinateOfSplittingList.size();i++)
	{
		printf("sorted splitting coordinate list: %f, with coordinate of splitting: %f with coordinate axis:%d\n",coordinateOfSplittingList[i], coordinateOfSplitting,nodeOfSubSpace->coordinateAxisOfSplitting);
	}
	*/
	nodeOfSubSpace->splittingPlaneCoordinate = coordinateOfSplitting;
	coordinateOfSplittingList.clear();
}

bool NearestTriangleSearch::terminateKDTreeConstruction(KDTreeNode * nodeOfSubSpace, unsigned int currentKDTreeDepth, unsigned int maxDepth, unsigned int minNumberOfTriangles){
	// abort: cell intersects with too few triangles and current tree is deep enough
	if ((nodeOfSubSpace->faceHandles->size()<=minNumberOfTriangles) || (currentKDTreeDepth==maxDepth)){ // TODO
		return true;
	}
	else{ 
		return false;
	}
}
}