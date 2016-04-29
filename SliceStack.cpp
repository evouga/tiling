#include "SliceStack.h"
#include "SliceParser.h"
#include "Slice.h"
#include "Tile.h"

#include <igl/viewer/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

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

void SliceStack::tetrahedralizeSlice (const Eigen::MatrixXd &botV, const Eigen::MatrixXi &botF, 
                                      const Eigen::MatrixXd &topV, const Eigen::MatrixXi &topF)
{
  Eigen::MatrixXd V(botV.rows() + topV.rows(), 3);
  Eigen::MatrixXi F(botF.rows() + topF.rows() + 8, 3);

  for (int i = 0; i < botV.rows(); ++i) {
    V.row(i) = botV.row(i);
  }

  for (int i = 0; i < topV.rows(); ++i) {
    V.row(i + botV.rows()) = topV.row(i);
  }

  int minXminYminZ = 0;
  int minXmaxYminZ = 3;
  int maxXminYminZ = 1;
  int maxXmaxYminZ = 2;

  int minXminYmaxZ = botV.rows() + minXminYminZ;
  int minXmaxYmaxZ = botV.rows() + minXmaxYminZ;
  int maxXminYmaxZ = botV.rows() + maxXminYminZ;
  int maxXmaxYmaxZ = botV.rows() + maxXmaxYminZ;

  for (int i = 0; i < botF.rows(); ++i) {
    F.row(i) = botF.row(i);
  }

  Eigen::Vector3i vertexOffset(botV.rows(), botV.rows(), botV.rows());

  for (int i = 0; i < topF.rows(); ++i) {
    F.row(i + botF.rows()) = vertexOffset + Eigen::Vector3i(topF.row(i));
  }

  int offset = botF.rows() + topF.rows();

  // left
  F.row(offset + 0) = Eigen::Vector3i(minXminYminZ, minXminYmaxZ, minXmaxYminZ);
  F.row(offset + 1) = Eigen::Vector3i(minXmaxYmaxZ, minXmaxYminZ, minXminYmaxZ);

  // right
  F.row(offset + 2) = Eigen::Vector3i(maxXminYminZ, maxXminYmaxZ, maxXmaxYminZ);
  F.row(offset + 3) = Eigen::Vector3i(maxXmaxYmaxZ, maxXmaxYminZ, maxXminYmaxZ);

  // front
  F.row(offset + 4) = Eigen::Vector3i(minXminYminZ, maxXminYminZ, minXminYmaxZ);
  F.row(offset + 5) = Eigen::Vector3i(maxXminYmaxZ, minXminYmaxZ, maxXminYminZ);

  // back
  F.row(offset + 6) = Eigen::Vector3i(maxXmaxYminZ, minXmaxYmaxZ, minXmaxYminZ);
  F.row(offset + 7) = Eigen::Vector3i(minXmaxYmaxZ, maxXmaxYminZ, maxXmaxYmaxZ);


  Eigen::MatrixXi testF(2, 3);
  testF.row(0) = Eigen::Vector3i(282,  262,  253);
  testF.row(1) = Eigen::Vector3i(156,    2,  155);

  cout << V.row(282) << endl << V.row(262) << endl << V.row(253) << endl;
  cout << V.row(156) << endl << V.row(2) << endl << V.row(155) << endl;
  
  cout << "bar " << endl;

  igl::writeOFF("foo.off", V, F);

  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, testF);
  viewer.data.set_face_based(true);
  viewer.launch();

  // Tetrahedralized interior
  // Eigen::MatrixXd TV;
  // Eigen::MatrixXi TT;
  // Eigen::MatrixXi TF;
  //
  // igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);

  cout << "Tetrahedralize done" << endl;
}
    
