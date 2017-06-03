/*
 *  MeshTraits.h
 *  BeNTO3D
 *
 *  Created by Hao Li on 2/15/08.
 *  Copyright 2008 ETH Zurich - Applied Geometry Group. All rights reserved.
 *
 */

#ifndef MESH_TRAITS_H
#define MESH_TRAITS_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

struct MeshTraits : public OpenMesh::DefaultTraits{

	typedef OpenMesh::Vec3f Point;
	typedef OpenMesh::Vec3f Normal;
	typedef OpenMesh::Vec2f TexCoord2D;
	typedef OpenMesh::Vec4f Color;
	
	VertexAttributes( OpenMesh::Attributes::Normal |
					  OpenMesh::Attributes::Color |
					  OpenMesh::Attributes::TexCoord2D);
					
	FaceAttributes( OpenMesh::Attributes::Normal |
                    OpenMesh::Attributes::Color );
						
	VertexTraits{
		private:
			bool selected;
			double doubleValue;
      int originalIndex;
		
		public:
			VertexT():selected(false),doubleValue(0.0),originalIndex(-1){}
		
			const bool& isSelected() const{return selected;}
			void setSelected(const bool _selected) {selected = _selected;}
			
			const double& getDoubleValue() const{return doubleValue;}
			void setDoubleValue(const double _doubleValue) {doubleValue = _doubleValue;}
			
			const double& getOriginalIndex() const{return originalIndex;}
			void setOriginalIndex(const double _originalIndex) {originalIndex = _originalIndex;}
		};
	
  EdgeTraits{
    private:
      bool _protected;
    public:
			EdgeT() : _protected(false){}
	
			const bool& isProtected() const{return _protected;}
			void setProtected(const bool _isProtected) {_protected = _isProtected;}
  };

	HalfedgeTraits{
  
		private:
			Color color;
			bool selected;

		public:
			HalfedgeT() : color( Color(0.0,0.0,0.0,1.0) ),selected(false){}

			const Color& get_color() const {return color;}
			void set_color(const Color& _color){color = _color;}

			const bool& isSelected() const{return selected;}
			void setSelected(const bool _selected) {selected = _selected;}
	};

	FaceTraits{
		private:
			bool selected;
		
		public:
			FaceT():selected(false){}
		
			const bool& isSelected() const{return selected;}
			void setSelected(const bool _selected) {selected = _selected;}
	};

};

#endif
