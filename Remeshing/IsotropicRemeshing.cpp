/*
 *  IsotropicRemeshing.cpp
 *  ScanRegistrator
 *
 *  Created by Hao Li on 02.06.07.
 *  Copyright 2007 ETH Zurich - Applied Geometry Group. All rights reserved.
 *
 */


#include "NearestTriangleSearch.h"
#include "IsotropicRemeshing.h"

#include <iostream>

using std::cout;
using std::endl;

IsotropicRemeshing::IsotropicRemeshing(double _targetLength)
{
	targetLength = _targetLength;
	lowerBoundLength = 4.0/5.0 * targetLength;
	higherBoundLength = 4.0/3.0 * targetLength;
}


IsotropicRemeshing::~IsotropicRemeshing()
{

}

void IsotropicRemeshing::remesh(TriangleMesh * _mesh,unsigned int _iterations,
                                ISOTROPIC_REMESHING_TYPES types)
{
	TriangleMesh mesh = *_mesh;
	removeNANVertices(&mesh);
	mesh.add_property(update);
	mesh.update_normals();
			
	meshBackup = mesh;
    for( unsigned int i = 0; i<_iterations; ++i )
	{
		splitLongEdges(&mesh);
		collapseShortEdges(&mesh);
		equalizeValences(&mesh);
		tangentialRelaxation(&mesh);	
		projectToSurface(&mesh);
	}
	areaEqualization(&mesh);
	mesh.remove_property(update);
	*_mesh = mesh;
}

bool IsotropicRemeshing::isNAN(double num)
{
	return num!=num;
}

void IsotropicRemeshing::removeNANVertices(TriangleMesh *_mesh)
{
	for(TriangleMesh::VIter vi = _mesh->vertices_begin(); vi != _mesh->vertices_end(); ++vi)
	{
		const TriangleMesh::Point &pt = _mesh->point( *vi );
		if(isNAN(pt[0]) || isNAN(pt[1]) || isNAN(pt[2]))
		{
			_mesh->delete_vertex( *vi, true );
		}
	}

	_mesh->garbage_collection();
}


void IsotropicRemeshing::splitLongEdges(TriangleMesh * _mesh)
{
	TriangleMesh::EIter     e_it, e_end;
	TriangleMesh::VHandle   v0, v1, vh;
	TriangleMesh::EHandle   eh, e0, e1;
	TriangleMesh::FHandle   f0, f1, f2, f3;
	bool            finished;
	int             i;



	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=_mesh->edges_begin(), e_end=_mesh->edges_end(); e_it!=e_end; ++e_it)
		{
      if (_mesh->data(*e_it).isProtected()) continue;

			v0 = _mesh->to_vertex_handle( _mesh->halfedge_handle( *e_it, 0 ) );
			v1 = _mesh->to_vertex_handle( _mesh->halfedge_handle( *e_it, 1 ) );

			const TriangleMesh::Point& p0 = _mesh->point(v0);
			const TriangleMesh::Point& p1 = _mesh->point(v1);

			if ( (p0-p1).norm() > higherBoundLength)
			{
				vh = _mesh->add_vertex((p0+p1)*0.5);
				_mesh->set_normal(vh, (_mesh->normal(v0) + _mesh->normal(v1)).normalize());
				_mesh->split( *e_it, vh );
				finished = false;
			}
		}
	}
}


void IsotropicRemeshing::collapseShortEdges(TriangleMesh * _mesh)
{
	TriangleMesh::EIter     e_it, e_end;
	TriangleMesh::CVVIter   vv_it;
	TriangleMesh::VHandle   v0, v1;
	TriangleMesh::HHandle   h0, h1, h01, h10;
	bool            finished, b0, b1;
	int             i;
	bool            hcol01, hcol10;

	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=_mesh->edges_begin(), e_end=_mesh->edges_end(); e_it!=e_end;++e_it)
		{
      if (_mesh->data(*e_it).isProtected()) continue;

            if( !_mesh->status( *e_it ).deleted() ) // might already be deleted
			{
				h10 = _mesh->halfedge_handle( *e_it, 0 );
				h01 = _mesh->halfedge_handle( *e_it, 1 );
				v0  = _mesh->to_vertex_handle(h10);
				v1  = _mesh->to_vertex_handle(h01);

				const TriangleMesh::Point& p0 = _mesh->point(v0);
				const TriangleMesh::Point& p1 = _mesh->point(v1);

				if ((p0-p1).length() < lowerBoundLength)
				{
					// get status
					b0 = _mesh->is_boundary(v0);
					b1 = _mesh->is_boundary(v1);
					hcol01 = hcol10 = true;


					// boundary rules (dont collapse boundary to interior)
					if( b0 && b1 ) { if( !_mesh->is_boundary( *e_it ) ) continue; }
					else if (b0) hcol01 = false; 
					else if (b1) hcol10 = false;
		
					// Testing before each collapse whether the collapse would produce an edge that is longer than high
					
					TriangleMesh::VertexVertexIter vvIt;
					vvIt=_mesh->vv_iter(v0);
					
					for(;vvIt.is_valid();++vvIt)
					{
						OpenMesh::Vec3f pointA = _mesh->point(v1);
						OpenMesh::Vec3f pointB = _mesh->point( *vvIt );
						if((pointA-pointB).length() >higherBoundLength)
						{
							hcol01 = false;
						}
					}		
					vvIt=_mesh->vv_iter(v1);
					for( ; vvIt.is_valid(); ++vvIt )
					{
						OpenMesh::Vec3f pointA = _mesh->point(v0);
						OpenMesh::Vec3f pointB = _mesh->point( *vvIt );
						if((pointA-pointB).length()>higherBoundLength)
						{
							hcol10 = false;
						}
					}													
		
					// topological rules
					if (hcol01)  hcol01 = _mesh->is_collapse_ok(h01);
					if (hcol10)  hcol10 = _mesh->is_collapse_ok(h10);


					// both collapses possible: collapse into vertex w/ higher valence
					
					if (hcol01 && hcol10)
					{
						if (_mesh->valence(v0) > _mesh->valence(v1))
							hcol10 = false;
						else
							hcol01 = false;
					}
					 

					// try v1 -> v0
					if (hcol10)
					{
						// contraint is too hard for boundaries
						//if(b1)
						//{
						//	_mesh->set_point(v0,_mesh->point(v1));
						//}
						_mesh->collapse(h10);
						finished = false;
					}

					// try v0 -> v1
					else if (hcol01)
					{
						//if(b0)
						//{
						//	_mesh->set_point(v1,_mesh->point(v0));
						//}
						_mesh->collapse(h01);
						finished = false;
					}
				}
			}
		}
	}

	_mesh->garbage_collection();

	if (i==100) std::cerr << "collapse break\n";
}

void IsotropicRemeshing::equalizeValences(TriangleMesh * _mesh)
{
  TriangleMesh::EIter     e_it, e_end;
  TriangleMesh::VHandle   v0, v1, v2, v3;
  TriangleMesh::HHandle   hh;
  int             val0, val1, val2, val3;
  int             val_opt0, val_opt1, val_opt2, val_opt3;
  int             ve0, ve1, ve2, ve3, ve_before, ve_after;
  bool            finished;
  int             i;
  
  // flip all edges
  for (finished=false, i=0; !finished && i<100; ++i)
  {
	  finished = true;

	  for (e_it=_mesh->edges_begin(), e_end=_mesh->edges_end(); e_it!=e_end; ++e_it)
	  {
		  if( !_mesh->is_boundary( *e_it ) && !_mesh->data(*e_it).isProtected())
		  {
			  hh = _mesh->halfedge_handle( *e_it, 0 );
			  v0 = _mesh->to_vertex_handle(hh);
			  v2 = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(hh));
			  //hh = _mesh->halfedge_handle(e_it, 1);
			  hh = _mesh->opposite_halfedge_handle(hh);
			  v1 = _mesh->to_vertex_handle(hh);
			  v3 = _mesh->to_vertex_handle(_mesh->next_halfedge_handle(hh));

			  val0 = _mesh->valence(v0);
			  val1 = _mesh->valence(v1);
			  val2 = _mesh->valence(v2);
			  val3 = _mesh->valence(v3);

			  val_opt0 = (_mesh->is_boundary(v0) ? 4 : 6);
			  val_opt1 = (_mesh->is_boundary(v1) ? 4 : 6);
			  val_opt2 = (_mesh->is_boundary(v2) ? 4 : 6);
			  val_opt3 = (_mesh->is_boundary(v3) ? 4 : 6);

			  ve0 = (val0 - val_opt0);
			  ve1 = (val1 - val_opt1);
			  ve2 = (val2 - val_opt2);
			  ve3 = (val3 - val_opt3);

			  ve_before = ve0*ve0 + ve1*ve1 + ve2*ve2 + ve3*ve3;
			  //ve_before = abs(ve0) + abs(ve1) + abs(ve2) + abs(ve3);
				
			  --val0;  --val1; 
			  ++val2;  ++val3;

			  ve0 = (val0 - val_opt0);
			  ve1 = (val1 - val_opt1);
			  ve2 = (val2 - val_opt2);
			  ve3 = (val3 - val_opt3);

			  ve_after = ve0*ve0 + ve1*ve1 + ve2*ve2 + ve3*ve3;
			  //ve_after = abs(ve0) + abs(ve1) + abs(ve2) + abs(ve3);

			  if( ve_before > ve_after && _mesh->is_flip_ok( *e_it ) )
			  {
				  _mesh->flip( *e_it );
				  finished = false;
			  }
		  }
	  }
  }

  if (i==100) std::cerr << "flip break\n";

}

void IsotropicRemeshing::tangentialRelaxation(TriangleMesh * _mesh)
{
	TriangleMesh::VIter     v_it, v_end(_mesh->vertices_end());
	TriangleMesh::CVVIter   vv_it;
	TriangleMesh::Scalar    valence;
	TriangleMesh::Point     u, n;

	
  // smooth
  for (int iters=0; iters<20; ++iters)
  {
    for (v_it=_mesh->vertices_begin(); v_it!=v_end; ++v_it)
    {
      if (!_mesh->is_boundary(*v_it))
      {
				u.vectorize(0.0);
				valence = 0;
				
				for( vv_it = _mesh->cvv_iter( *v_it ); vv_it.is_valid(); ++vv_it )
				{
					u += _mesh->point( *vv_it );
					++valence;
				}
				
				if (valence)
				{
					u *= (float)(1.0/valence);
					u -= _mesh->point( *v_it );
					n = _mesh->normal( *v_it );
					n *= (u | n);
					u -= n;
				}
				
				_mesh->property( update, *v_it ) = u;
      }
    }
		
    for (v_it=_mesh->vertices_begin(); v_it!=v_end; ++v_it)
		if( !_mesh->is_boundary( *v_it ) )
			_mesh->point( *v_it ) += _mesh->property( update, *v_it );
  }
}

void IsotropicRemeshing::projectToSurface(TriangleMesh * _mesh)
{

	dataStructure::NearestTriangleSearch * triangleKDTree = new dataStructure::NearestTriangleSearch();
	triangleKDTree->initializeKDTreeBasedSearchStructure(&meshBackup,12,20);
	
	TriangleMesh::VertexIter vIt;
		
	// iterate over vertices of mesh
	for(vIt=_mesh->vertices_begin();vIt!=_mesh->vertices_end();++vIt)
	{
		if( !_mesh->is_boundary( *vIt ) )
		{
			TriangleMesh::FaceHandle closestTriangleOnBackupMesh;
			OpenMesh::Vec3f pointOnMesh;
			pointOnMesh = _mesh->point( *vIt );
			assert(!isNAN(pointOnMesh[0]));
			assert(!isNAN(pointOnMesh[1]));
			assert(!isNAN(pointOnMesh[2]));
			triangleKDTree->computeClosestTriangleOfPoint(pointOnMesh,&meshBackup,closestTriangleOnBackupMesh);
			// project point of mesh to closest triangle on backup mesh
			OpenMesh::Vec3f pointOnBackupMeshFromMesh;
			TriangleMesh::FaceVertexIter fvIt;
			
			fvIt = meshBackup.fv_iter(closestTriangleOnBackupMesh);

			// Shouldn't be needed, but included just to be safe...
			if(!closestTriangleOnBackupMesh.is_valid())
				continue;
			
			OpenMesh::Vec3f pointA = meshBackup.point( *fvIt );
			OpenMesh::Vec3f pointB = meshBackup.point( *(++fvIt) );
			OpenMesh::Vec3f pointC = meshBackup.point( *(++fvIt) );
			
			
			geometry::squaredDistancePointToTriangle(	pointOnMesh,
														pointOnBackupMeshFromMesh,
														pointA,
														pointB,
														pointC);
																																		
			_mesh->set_point( *vIt, pointOnBackupMeshFromMesh );
		}
	}
	
	
	_mesh->update_face_normals();
	_mesh->update_vertex_normals();
	
	triangleKDTree->destroyKDTreeBasedSearchStructure();
}


void IsotropicRemeshing::areaEqualization(TriangleMesh * _mesh)
{
	TriangleMesh::VIter     v_it, v_end(_mesh->vertices_end());
	TriangleMesh::CVVIter   vv_it;
	TriangleMesh::Scalar    valence;
	TriangleMesh::Point     u, n;
	
	for (int iters=0; iters<5; ++iters)
	{
		// smooth
		for (v_it=_mesh->vertices_begin(); v_it!=v_end; ++v_it)
		{
			if( !_mesh->is_boundary( *v_it ) )
		  {
					u.vectorize(0.0);
					valence = 0;
					
					double totalNeighborAreas= 0.0;
					for( vv_it = _mesh->cvv_iter( *v_it ); vv_it.is_valid(); ++vv_it )
					{
						TriangleMesh::VFIter vfIt;
						// compute neighbor area
						
						double neighborArea = 0.0;
						for( vfIt = _mesh->vf_iter( *vv_it ); vfIt.is_valid(); ++vfIt )
						{
							OpenMesh::FaceHandle fh= *vfIt;
							neighborArea += geometry::triangleSurfaceArea(fh,_mesh);
						}
						totalNeighborAreas += neighborArea;
						u += _mesh->point(*vv_it)*neighborArea;
						++valence;
					}
					
					
					if (valence && totalNeighborAreas != 0.0)
					{
						double lambda = 0.5;
						u *= (float)(1.0/totalNeighborAreas);
						u -= _mesh->point( *v_it );
						n = _mesh->normal( *v_it );
						n *= (u | n);
						u -= n*lambda;
					}
					
					_mesh->property( update, *v_it ) = u;
		  }
		}
			
		for (v_it=_mesh->vertices_begin(); v_it!=v_end; ++v_it)
			if( !_mesh->is_boundary( *v_it ) )
				_mesh->point( *v_it ) += _mesh->property( update, *v_it );
	  
	}
}
