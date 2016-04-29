#include <vector>
#include <algorithm>

#include "SliceStack.h"
#include "SliceParser.h"
#include "Slice.h"
#include "Tile.h"

#include <igl/viewer/Viewer.h>
#include <igl/triangle/triangulate.h>
// #include <igl/copyleft/tetgen/tetrahedralize.h>
// #include <igl/copyleft/cgal/remesh_self_intersections.h>

using namespace std;

SliceStack::SliceStack(const char *baseFilename, const char *objectname)
{
  readSlicesFromFolder(baseFilename, objectname, slices_);
  numSlices_ = slices_.size();
  double z = 0;
  for(int i=0; i<numSlices_; i++)
  {
    heights_.push_back(z);
    z += slices_[i]->thickness;
  }
}

SliceStack::~SliceStack()
{
  for(int i=0; i<numSlices_; i++)
    delete slices_[i];
}

void SliceStack::triangulateSlice(int bottomidx, double areaBound, 
                                  Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
                                  Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces)
{
  assert(bottomidx >= 0 && bottomidx < slices_.size());
  Tile t(*slices_[bottomidx], *slices_[bottomidx+1]);
  t.triangulateSlices(areaBound, botverts, botfaces, topverts, topfaces);
}

int SliceStack::getSizeAt(int i) { 
  if (i > slices_.size()) return -1;
  return slices_[i]->contours.size();
}

bool customSortByX (const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[0] < b[0];
}

bool customSortByY (const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[1] < b[1];
}

bool customSortByZ (const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[2] < b[2];
}

// constantCoord => 0: (-0.5, _, _), 1: (0.5, _, _), 2: (_, -0.5, _), 3: (_, 0.5, _)
void SliceStack::triangulateSide (int constantCoord, vector<Eigen::Vector3d> &verts, 
                                  Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
  sort(verts.begin(), verts.end(), customSortByZ);

  vector<Eigen::Vector3d> topV;
  vector<Eigen::Vector3d> botV;

  for (int i = 0; i < verts.size(); ++i) {
    if (verts[i][2] < 0.0)
      botV.push_back(verts[i]);
    else
      topV.push_back(verts[i]);
  }

  if (constantCoord >= 2) {
    sort(botV.begin(), botV.end(), customSortByX);
    sort(topV.rbegin(), topV.rend(), customSortByX);
  }
  else {
    sort(botV.begin(), botV.end(), customSortByY);
    sort(topV.rbegin(), topV.rend(), customSortByY);
  }

  Eigen::MatrixXd inputV(verts.size(), 2);
  Eigen::MatrixXi inputE(verts.size(), 2);
  Eigen::MatrixXd H(0, 2);

  for (int i = 0; i < botV.size(); ++i) {
    if (constantCoord < 2)
      inputV.row(i) = Eigen::Vector2d(botV[i][1], botV[i][2]);
    else 
      inputV.row(i) = Eigen::Vector2d(botV[i][0], botV[i][2]);
  }

  for (int i = 0; i < topV.size(); ++i) {
    if (constantCoord >= 2)
      inputV(botV.size() + i, 0) = topV[i][0];
    else
      inputV(botV.size() + i, 0) = topV[i][1];
    inputV(botV.size() + i, 1) = topV[i][2];
  }

  for (int i = 0; i < verts.size(); ++i) {
    inputE.row(i) = Eigen::Vector2i(i, (i + 1) % (verts.size()));
  }

  Eigen::MatrixXd tmpV;

	cout << "Delaunay triangulating mesh with " << inputV.rows() << " verts" << endl;
	igl::triangle::triangulate(inputV,inputE,H,"Da0.05q",tmpV,F);	
	cout << "done" << endl;

  V.resize(tmpV.rows(), 3);

  for (int i = 0; i < tmpV.rows(); ++i) {
    if (constantCoord == 0)
      V.row(i) = Eigen::Vector3d(-0.5, tmpV(i, 0), tmpV(i, 1));
    else if (constantCoord == 1)
      V.row(i) = Eigen::Vector3d(0.5, tmpV(i, 0), tmpV(i, 1));
    else if (constantCoord == 2)
      V.row(i) = Eigen::Vector3d(tmpV(i, 0), -0.5, tmpV(i, 1));
    else if (constantCoord == 3)
      V.row(i) = Eigen::Vector3d(tmpV(i, 0), 0.5, tmpV(i, 1));
  }
}

void SliceStack::tetrahedralizeSlice (const Eigen::MatrixXd &botV, const Eigen::MatrixXi &botF, 
                                      const Eigen::MatrixXd &topV, const Eigen::MatrixXi &topF)
{
  vector<Eigen::Vector3d> leftV;
  vector<Eigen::Vector3d> rightV;
  vector<Eigen::Vector3d> frontV;
  vector<Eigen::Vector3d> backV;

  double bound = 0.5;

  for (int i = 0; i < botV.rows(); ++i) {
    // x coordinate touches
    if (botV(i, 0) == -bound) {
      leftV.push_back(botV.row(i));
    }
    if (botV(i, 0) == bound) {
      rightV.push_back(botV.row(i));
    }
    // y coordinate touches
    if (botV(i, 1) == -bound) {
      frontV.push_back(botV.row(i)); 
    }
    if (botV(i, 1) == bound) {
      backV.push_back(botV.row(i)); 
    }
  }

  for (int i = 0; i < topV.rows(); ++i) {
    // x coordinate touches
    if (topV(i, 0) == -bound) {
      leftV.push_back(topV.row(i));
    }
    if (topV(i, 0) == bound) {
      rightV.push_back(topV.row(i));
    }
    // y coordinate touches
    if (topV(i, 1) == -bound) {
      frontV.push_back(topV.row(i)); 
    }
    if (topV(i, 1) == bound) {
      backV.push_back(topV.row(i)); 
    }
  }

  Eigen::MatrixXd leftTriV;
  Eigen::MatrixXi leftTriF;

  Eigen::MatrixXd rightTriV;
  Eigen::MatrixXi rightTriF;

  Eigen::MatrixXd backTriV;
  Eigen::MatrixXi backTriF;

  Eigen::MatrixXd frontTriV;
  Eigen::MatrixXi frontTriF;

  triangulateSide(0, leftV, leftTriV, leftTriF);

  triangulateSide(2, frontV, frontTriV, frontTriF);

  triangulateSide(1, rightV, rightTriV, rightTriF);

  triangulateSide(3, backV, backTriV, backTriF);

  int totalVertices = topV.rows() + botV.rows() + 
    leftTriV.rows() + rightTriV.rows() + backTriV.rows() + frontTriV.rows();

  int totalFaces = topF.rows() + botF.rows() + 
    leftTriF.rows() + rightTriF.rows() + backTriF.rows() + frontTriF.rows();

  Eigen::MatrixXd V(totalVertices, 3);
  Eigen::MatrixXi F(totalFaces, 3);

  int offset = 0;

  // Add the vertices
  for (int i = 0; i < botV.rows(); ++i) {
    V.row(i) = botV.row(i);
  }

  offset += botV.rows();

  for (int i = 0; i < topV.rows(); ++i) {
    V.row(offset + i) = topV.row(i);
  }

  offset += topV.rows();

  for (int i = 0; i < leftTriV.rows(); ++i) {
    V.row(offset + i) = leftTriV.row(i);
  }

  offset += leftTriV.rows();

  for (int i = 0; i < rightTriV.rows(); ++i) {
    V.row(offset + i) = rightTriV.row(i);
  }

  offset += rightTriV.rows();

  for (int i = 0; i < frontTriV.rows(); ++i) {
    V.row(offset + i) = frontTriV.row(i);
  }

  offset += frontTriV.rows();

  for (int i = 0; i < backTriV.rows(); ++i) {
    V.row(offset + i) = backTriV.row(i);
  }

  // Add the faces
  Eigen::Vector3i vertexOffset(0, 0, 0);
  offset = 0;

  for (int i = 0; i < botF.rows(); ++i) {
    F.row(i) = botF.row(i);
  }

  vertexOffset += Eigen::Vector3i(botV.rows(), botV.rows(), botV.rows());
  offset += botF.rows();

  for (int i = 0; i < topF.rows(); ++i) {
    F.row(offset + i) = vertexOffset + Eigen::Vector3i(topF.row(i));
  }

  vertexOffset += Eigen::Vector3i(topV.rows(), topV.rows(), topV.rows());
  offset += topF.rows();

  for (int i = 0; i < leftTriF.rows(); ++i) {
    F.row(offset + i) = vertexOffset + Eigen::Vector3i(leftTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(leftTriV.rows(), leftTriV.rows(), leftTriV.rows());
  offset += leftTriF.rows();

  for (int i = 0; i < rightTriF.rows(); ++i) {
    F.row(offset + i) = vertexOffset + Eigen::Vector3i(rightTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(rightTriV.rows(), rightTriV.rows(), rightTriV.rows());
  offset += rightTriF.rows();

  for (int i = 0; i < frontTriF.rows(); ++i) {
    F.row(offset + i) = vertexOffset + Eigen::Vector3i(frontTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(frontTriV.rows(), frontTriV.rows(), frontTriV.rows());
  offset += frontTriF.rows();

  for (int i = 0; i < backTriF.rows(); ++i) {
    F.row(offset + i) = vertexOffset + Eigen::Vector3i(backTriF.row(i));
  }

  cout << "bar " << endl;
  
  // igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
  // Eigen::MatrixXd VV;
  // Eigen::MatrixXi FF;
  // Eigen::MatrixXi IF;
  // Eigen::VectorXi J, IM;
  //
  // igl::copyleft::cgal::remesh_self_intersections(V,F,params,VV,FF,IF,J,IM);
  //
  // igl::writeOFF("foo.off", VV, FF);

  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_colors(colors);
  viewer.data.set_face_based(false);
  viewer.launch();

  // Tetrahedralized interior
  // Eigen::MatrixXd TV;
  // Eigen::MatrixXi TT;
  // Eigen::MatrixXi TF;
  //
  // igl::copyleft::tetgen::tetrahedralize(VV,FF,"pqd1.414Y", TV,TT,TF);

  cout << "Tetrahedralize done" << endl;
}
