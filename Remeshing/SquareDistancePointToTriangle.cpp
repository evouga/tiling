#include "BasicMeasurements.h"

namespace geometry{

double squaredDistancePointToTriangle(	TriangleMesh::Point & _point,
										TriangleMesh::Point & _triangleVertexA,
										TriangleMesh::Point & _triangleVertexB,
										TriangleMesh::Point & _triangleVertexC)
{

  TriangleMesh::Point kDiff = _triangleVertexA - _point;
  TriangleMesh::Point kEdge0 = _triangleVertexB - _triangleVertexA;
  TriangleMesh::Point kEdge1 = _triangleVertexC - _triangleVertexA;
 double fA00 = kEdge0.sqrnorm();
 double fA01 = kEdge0 | kEdge1;
 double fA11 = kEdge1.sqrnorm();
 double fB0 = kDiff | kEdge0;
 double fB1 = kDiff | kEdge1;
 double fC = kDiff.sqrnorm();
 double fDet = fabs(fA00*fA11-fA01*fA01);
 double fS = fA01*fB1-fA11*fB0;
 double fT = fA01*fB0-fA00*fB1;
 double fSqrDistance;

 if (fS + fT <= fDet)
 {
     if (fS < (double)0.0)
     {
         if (fT < (double)0.0)  // region 4
         {
             if (fB0 < (double)0.0)
             {
                 fT = (double)0.0;
                 if (-fB0 >= fA00)
                 {
                     fS = (double)1.0;
                     fSqrDistance = fA00+((double)2.0)*fB0+fC;
                 }
                 else
                 {
                     fS = -fB0/fA00;
                     fSqrDistance = fB0*fS+fC;
                 }
             }
             else
             {
                 fS = (double)0.0;
                 if (fB1 >= (double)0.0)
                 {
                     fT = (double)0.0;
                     fSqrDistance = fC;
                 }
                 else if (-fB1 >= fA11)
                 {
                     fT = (double)1.0;
                     fSqrDistance = fA11+((double)2.0)*fB1+fC;
                 }
                 else
                 {
                     fT = -fB1/fA11;
                     fSqrDistance = fB1*fT+fC;
                 }
             }
         }
         else  // region 3
         {
             fS = (double)0.0;
             if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else if (-fB1 >= fA11)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 5
     {
         fT = (double)0.0;
         if (fB0 >= (double)0.0)
         {
             fS = (double)0.0;
             fSqrDistance = fC;
         }
         else if (-fB0 >= fA00)
         {
             fS = (double)1.0;
             fSqrDistance = fA00+((double)2.0)*fB0+fC;
         }
         else
         {
             fS = -fB0/fA00;
             fSqrDistance = fB0*fS+fC;
         }
     }
     else  // region 0
     {
         // minimum at interior point
         double fInvDet = ((double)1.0)/fDet;
         fS *= fInvDet;
         fT *= fInvDet;
         fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
             fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
     }
 }
 else
 {
     double fTmp0, fTmp1, fNumer, fDenom;

     if (fS < (double)0.0)  // region 2
     {
         fTmp0 = fA01 + fB0;
         fTmp1 = fA11 + fB1;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fS = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 6
     {
         fTmp0 = fA01 + fB1;
         fTmp1 = fA00 + fB0;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-((double)2.0)*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fT = (double)1.0;
                 fS = (double)0.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = fNumer/fDenom;
                 fS = (double)1.0 - fT;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fT = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fS = (double)1.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else if (fB0 >= (double)0.0)
             {
                 fS = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fS = -fB0/fA00;
                 fSqrDistance = fB0*fS+fC;
             }
         }
     }
     else  // region 1
     {
         fNumer = fA11 + fB1 - fA01 - fB0;
         if (fNumer <= (double)0.0)
         {
             fS = (double)0.0;
             fT = (double)1.0;
             fSqrDistance = fA11+((double)2.0)*fB1+fC;
         }
         else
         {
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
     }
  }

  // account for numerical round-off error
  if (fSqrDistance < (double)0.0)
  {
     fSqrDistance = (double)0.0;
  }

  return fSqrDistance; 
}

double squaredDistancePointToTriangle(			TriangleMesh::Point & _point,
												TriangleMesh::Point & _projectedPoint,
												TriangleMesh::Point & _triangleVertexA,
												TriangleMesh::Point & _triangleVertexB,
												TriangleMesh::Point & _triangleVertexC)
{
  TriangleMesh::Point kDiff = _triangleVertexA - _point;
  TriangleMesh::Point kEdge0 = _triangleVertexB - _triangleVertexA;
  TriangleMesh::Point kEdge1 = _triangleVertexC - _triangleVertexA;
 double fA00 = kEdge0.sqrnorm();
 double fA01 = kEdge0 | kEdge1;
 double fA11 = kEdge1.sqrnorm();
 double fB0 = kDiff | kEdge0;
 double fB1 = kDiff | kEdge1;
 double fC = kDiff.sqrnorm();
 double fDet = fabs(fA00*fA11-fA01*fA01);
 double fS = fA01*fB1-fA11*fB0;
 double fT = fA01*fB0-fA00*fB1;
 double fSqrDistance;

 if (fS + fT <= fDet)
 {
     if (fS < (double)0.0)
     {
         if (fT < (double)0.0)  // region 4
         {
             if (fB0 < (double)0.0)
             {
                 fT = (double)0.0;
                 if (-fB0 >= fA00)
                 {
                     fS = (double)1.0;
                     fSqrDistance = fA00+((double)2.0)*fB0+fC;
                 }
                 else
                 {
                     fS = -fB0/fA00;
                     fSqrDistance = fB0*fS+fC;
                 }
             }
             else
             {
                 fS = (double)0.0;
                 if (fB1 >= (double)0.0)
                 {
                     fT = (double)0.0;
                     fSqrDistance = fC;
                 }
                 else if (-fB1 >= fA11)
                 {
                     fT = (double)1.0;
                     fSqrDistance = fA11+((double)2.0)*fB1+fC;
                 }
                 else
                 {
                     fT = -fB1/fA11;
                     fSqrDistance = fB1*fT+fC;
                 }
             }
         }
         else  // region 3
         {
             fS = (double)0.0;
             if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else if (-fB1 >= fA11)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 5
     {
         fT = (double)0.0;
         if (fB0 >= (double)0.0)
         {
             fS = (double)0.0;
             fSqrDistance = fC;
         }
         else if (-fB0 >= fA00)
         {
             fS = (double)1.0;
             fSqrDistance = fA00+((double)2.0)*fB0+fC;
         }
         else
         {
             fS = -fB0/fA00;
             fSqrDistance = fB0*fS+fC;
         }
     }
     else  // region 0
     {
         // minimum at interior point
         double fInvDet = ((double)1.0)/fDet;
         fS *= fInvDet;
         fT *= fInvDet;
         fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
             fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
     }
 }
 else
 {
     double fTmp0, fTmp1, fNumer, fDenom;

     if (fS < (double)0.0)  // region 2
     {
         fTmp0 = fA01 + fB0;
         fTmp1 = fA11 + fB1;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fS = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 6
     {
         fTmp0 = fA01 + fB1;
         fTmp1 = fA00 + fB0;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-((double)2.0)*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fT = (double)1.0;
                 fS = (double)0.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = fNumer/fDenom;
                 fS = (double)1.0 - fT;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fT = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fS = (double)1.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else if (fB0 >= (double)0.0)
             {
                 fS = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fS = -fB0/fA00;
                 fSqrDistance = fB0*fS+fC;
             }
         }
     }
     else  // region 1
     {
         fNumer = fA11 + fB1 - fA01 - fB0;
         if (fNumer <= (double)0.0)
         {
             fS = (double)0.0;
             fT = (double)1.0;
             fSqrDistance = fA11+((double)2.0)*fB1+fC;
         }
         else
         {
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
     }
  }

  // account for numerical round-off error
  if (fSqrDistance < (double)0.0)
  {
     fSqrDistance = (double)0.0;
  }	
	
	_projectedPoint = _triangleVertexA + (TriangleMesh::Scalar)fS*kEdge0 + (TriangleMesh::Scalar)fT*kEdge1;
	return fSqrDistance;
}


double squaredDistancePointToTriangle(			TriangleMesh::Point & _point,
												TriangleMesh::Point & _projectedPoint,
												TriangleMesh::Point & _projectedPointBaryCentricCoordinates,
												TriangleMesh::Point & _triangleVertexA,
												TriangleMesh::Point & _triangleVertexB,
												TriangleMesh::Point & _triangleVertexC
												)
{
  TriangleMesh::Point kDiff = _triangleVertexA - _point;
  TriangleMesh::Point kEdge0 = _triangleVertexB - _triangleVertexA;
  TriangleMesh::Point kEdge1 = _triangleVertexC - _triangleVertexA;
 double fA00 = kEdge0.sqrnorm();
 double fA01 = kEdge0 | kEdge1;
 double fA11 = kEdge1.sqrnorm();
 double fB0 = kDiff | kEdge0;
 double fB1 = kDiff | kEdge1;
 double fC = kDiff.sqrnorm();
 double fDet = fabs(fA00*fA11-fA01*fA01);
 double fS = fA01*fB1-fA11*fB0;
 double fT = fA01*fB0-fA00*fB1;
 double fSqrDistance;

 if (fS + fT <= fDet)
 {
     if (fS < (double)0.0)
     {
         if (fT < (double)0.0)  // region 4
         {
             if (fB0 < (double)0.0)
             {
                 fT = (double)0.0;
                 if (-fB0 >= fA00)
                 {
                     fS = (double)1.0;
                     fSqrDistance = fA00+((double)2.0)*fB0+fC;
                 }
                 else
                 {
                     fS = -fB0/fA00;
                     fSqrDistance = fB0*fS+fC;
                 }
             }
             else
             {
                 fS = (double)0.0;
                 if (fB1 >= (double)0.0)
                 {
                     fT = (double)0.0;
                     fSqrDistance = fC;
                 }
                 else if (-fB1 >= fA11)
                 {
                     fT = (double)1.0;
                     fSqrDistance = fA11+((double)2.0)*fB1+fC;
                 }
                 else
                 {
                     fT = -fB1/fA11;
                     fSqrDistance = fB1*fT+fC;
                 }
             }
         }
         else  // region 3
         {
             fS = (double)0.0;
             if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else if (-fB1 >= fA11)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 5
     {
         fT = (double)0.0;
         if (fB0 >= (double)0.0)
         {
             fS = (double)0.0;
             fSqrDistance = fC;
         }
         else if (-fB0 >= fA00)
         {
             fS = (double)1.0;
             fSqrDistance = fA00+((double)2.0)*fB0+fC;
         }
         else
         {
             fS = -fB0/fA00;
             fSqrDistance = fB0*fS+fC;
         }
     }
     else  // region 0
     {
         // minimum at interior point
         double fInvDet = ((double)1.0)/fDet;
         fS *= fInvDet;
         fT *= fInvDet;
         fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
             fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
     }
 }
 else
 {
     double fTmp0, fTmp1, fNumer, fDenom;

     if (fS < (double)0.0)  // region 2
     {
         fTmp0 = fA01 + fB0;
         fTmp1 = fA11 + fB1;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+2.0f*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fS = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fT = (double)1.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else if (fB1 >= (double)0.0)
             {
                 fT = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fT = -fB1/fA11;
                 fSqrDistance = fB1*fT+fC;
             }
         }
     }
     else if (fT < (double)0.0)  // region 6
     {
         fTmp0 = fA01 + fB1;
         fTmp1 = fA00 + fB0;
         if (fTmp1 > fTmp0)
         {
             fNumer = fTmp1 - fTmp0;
             fDenom = fA00-((double)2.0)*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fT = (double)1.0;
                 fS = (double)0.0;
                 fSqrDistance = fA11+((double)2.0)*fB1+fC;
             }
             else
             {
                 fT = fNumer/fDenom;
                 fS = (double)1.0 - fT;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
         else
         {
             fT = (double)0.0;
             if (fTmp1 <= (double)0.0)
             {
                 fS = (double)1.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else if (fB0 >= (double)0.0)
             {
                 fS = (double)0.0;
                 fSqrDistance = fC;
             }
             else
             {
                 fS = -fB0/fA00;
                 fSqrDistance = fB0*fS+fC;
             }
         }
     }
     else  // region 1
     {
         fNumer = fA11 + fB1 - fA01 - fB0;
         if (fNumer <= (double)0.0)
         {
             fS = (double)0.0;
             fT = (double)1.0;
             fSqrDistance = fA11+((double)2.0)*fB1+fC;
         }
         else
         {
             fDenom = fA00-2.0f*fA01+fA11;
             if (fNumer >= fDenom)
             {
                 fS = (double)1.0;
                 fT = (double)0.0;
                 fSqrDistance = fA00+((double)2.0)*fB0+fC;
             }
             else
             {
                 fS = fNumer/fDenom;
                 fT = (double)1.0 - fS;
                 fSqrDistance = fS*(fA00*fS+fA01*fT+((double)2.0)*fB0) +
                     fT*(fA01*fS+fA11*fT+((double)2.0)*fB1)+fC;
             }
         }
     }
  }

  // account for numerical round-off error
  if (fSqrDistance < (double)0.0)
  {
     fSqrDistance = (double)0.0;
  }	
	
	_projectedPoint = _triangleVertexA + (TriangleMesh::Scalar)fS*kEdge0 + (TriangleMesh::Scalar)fT*kEdge1;
	
	_projectedPointBaryCentricCoordinates[0]=1.0-fS-fT;
	_projectedPointBaryCentricCoordinates[1]=fS;
	_projectedPointBaryCentricCoordinates[2]=fT;	
	
	return fSqrDistance;
}




}



