#include "BasicMeasurements.h"

namespace geometry{

double xsqrt(double value)
{
	return value > 0. ? sqrt(value) : 0.;
}

double triangleSurfaceArea(	TriangleMesh::FaceHandle & _face,
							TriangleMesh * inputMesh)
{
	
	TriangleMesh::FaceVertexIter fvIt=inputMesh->fv_iter(_face);
	OpenMesh::Vec3f pointA = inputMesh->point( *fvIt );
	OpenMesh::Vec3f pointB = inputMesh->point( *(++fvIt) );
	OpenMesh::Vec3f pointC = inputMesh->point( *(++fvIt) );

	double a = (pointA-pointB).length();
	double b = (pointB-pointC).length();
	double c = (pointC-pointA).length();

	if(a>=b){
		if(b>=c){
			// a>=b>=c
			return 0.25*xsqrt(
				(a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
				);
		}
		else{
			if(c>=a){
				// c>=a>=b
				return 0.25*xsqrt(
				(c+(a+b))*(b-(c-a))*(b+(c-a))*(c+(a-b))
				);
			}
			else{
				// a>=c>=b
				return 0.25*xsqrt(
				(a+(c+b))*(b-(a-c))*(b+(a-c))*(a+(c-b))
				);

			}
		}
	}
	else{
		if(a>=c){
			// b>a>=c
			return 0.25*xsqrt(
				(b+(c+a))*(a-(b-c))*(a+(b-c))*(b+(c-a))
				);
		}
		else{
			if(c>=b){
				// c>=b>=a
				return 0.25*xsqrt(
				(c+(b+a))*(a-(c-b))*(a+(c-b))*(c+(b-a))
				);
			}
			else{
				// b>=c>=a
				return 0.25*xsqrt(
				(b+(a+c))*(c-(b-a))*(c+(b-a))*(b+(a-c))
				);

			}
		}
	}
}

double triangleSurfaceArea(	OpenMesh::Vec3f & pointA,
							OpenMesh::Vec3f & pointB,
							OpenMesh::Vec3f & pointC){

	double a = (pointA-pointB).length();
	double b = (pointB-pointC).length();
	double c = (pointC-pointA).length();

	if(a>=b){
		if(b>=c){
			// a>=b>=c
			return 0.25*xsqrt(
				(a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c))
				);
		}
		else{
			if(c>=a){
				// c>=a>=b
				return 0.25*xsqrt(
				(c+(a+b))*(b-(c-a))*(b+(c-a))*(c+(a-b))
				);
			}
			else{
				// a>=c>=b
				return 0.25*xsqrt(
				(a+(c+b))*(b-(a-c))*(b+(a-c))*(a+(c-b))
				);

			}
		}
	}
	else{
		if(a>=c){
			// b>a>=c
			return 0.25*xsqrt(
				(b+(c+a))*(a-(b-c))*(a+(b-c))*(b+(c-a))
				);
		}
		else{
			if(c>=b){
				// c>=b>=a
				return 0.25*xsqrt(
				(c+(b+a))*(a-(c-b))*(a+(c-b))*(c+(b-a))
				);
			}
			else{
				// b>=c>=a
				return 0.25*xsqrt(
				(b+(a+c))*(c-(b-a))*(c+(b-a))*(b+(a-c))
				);

			}
		}
	}
}

}