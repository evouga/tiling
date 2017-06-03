/*
 *  IsotropicRemeshing.h
 *  ScanRegistrator
 *
 *  Created by Hao Li on 02.06.07.
 *  Copyright 2007 ETH Zurich - Applied Geometry Group. All rights reserved.
 *
 */

#ifndef ISOTROPIC_REMESHING_H
#define ISOTROPIC_REMESHING_H

#include "MeshTraits.h"

typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> TriangleMesh;
enum ISOTROPIC_REMESHING_TYPES {
  SPLIT_LONG_EDGES = 0x0001,
  COLLAPSE_SHORT_EDGES = 0x0002,
  EQUALIZE_VALENCES = 0x0004,
  TANGENTIAL_RELAXATION = 0x0008,
  PROJECT_TO_SURFACE = 0x0010,
  ALL = 0x001f
};

class IsotropicRemeshing
{
	public:
	
    IsotropicRemeshing(double _targetLength);
	~IsotropicRemeshing();
	
	void remesh(TriangleMesh * _mesh,unsigned int _iterations, 
              ISOTROPIC_REMESHING_TYPES types = ALL);
	
	protected:
		
	bool isNAN(double num);
	void removeNANVertices(TriangleMesh * _mesh);
	void splitLongEdges(TriangleMesh * _mesh);	
	void collapseShortEdges(TriangleMesh * _mesh);
	void equalizeValences(TriangleMesh * _mesh);
	void tangentialRelaxation(TriangleMesh * _mesh);
	void projectToSurface(TriangleMesh * _mesh);
	
	void areaEqualization(TriangleMesh * _mesh);
	
	double targetLength;
	double lowerBoundLength;
	double higherBoundLength;
	
	TriangleMesh meshBackup;
	
	OpenMesh::VPropHandleT<TriangleMesh::Point> update;
};

#endif //ISOTROPIC_REMESHING_H 
