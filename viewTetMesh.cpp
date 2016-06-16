#include "viewTetMesh.h"

#include <igl/barycenter.h>
#include <igl/bfs_orient.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/embree/reorient_facets_raycast.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>

// These need to be global, so they can be used in Viewer functions
namespace {
  // Barycenters, used for defining tet centers
  Eigen::MatrixXd _B; // Internal barycenters
  Eigen::MatrixXd _C; // Internal color values
  // Tetrahedralized interior, passed as parameters
  Eigen::MatrixXd _TV; // tet vertices
  Eigen::MatrixXi _TT; // tet tets (similar to triangles)
  Eigen::MatrixXi _TF; // tet faces
  Eigen::VectorXd _TC; // tet vertex values (will be cerated to colors)
} // namespace


namespace {
// This function is called every time a keyboard button is pressed
bool key_down_depth(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key >= '1' && key <= '9')
  {
    double t = double((key - '1')+1) / 9.0;
		double t_abs = _TV.col(2).minCoeff() + 
				(_TV.col(2).maxCoeff() - _TV.col(2).minCoeff()) * t;

    VectorXd v = _B.col(2).array() - _B.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    vector<int> s;

    for (unsigned i=0; i<v.size(); ++i) {
      if (v(i) <= t) {
        s.push_back(i);
			}
		}

    MatrixXd V_temp(s.size()*4,3);
    MatrixXi F_temp(s.size()*4,3);
		MatrixXd C_temp(s.size()*4,3);
    VectorXd S_temp(s.size()*4);

    for (unsigned i=0; i<s.size();++i)
    {
      V_temp.row(i*4+0) = _TV.row(_TT(s[i],0));
      V_temp.row(i*4+1) = _TV.row(_TT(s[i],1));
      V_temp.row(i*4+2) = _TV.row(_TT(s[i],2));
      V_temp.row(i*4+3) = _TV.row(_TT(s[i],3));
			if (_TC.rows() > 0) {
				C_temp.row(i*4+0) = _C.row(_TT(s[i],0));
				C_temp.row(i*4+1) = _C.row(_TT(s[i],1));
				C_temp.row(i*4+2) = _C.row(_TT(s[i],2));
				C_temp.row(i*4+3) = _C.row(_TT(s[i],3));
        S_temp(i*4+0) = _TC(_TT(s[i],0));
        S_temp(i*4+1) = _TC(_TT(s[i],1));
        S_temp(i*4+2) = _TC(_TT(s[i],2));
        S_temp(i*4+3) = _TC(_TT(s[i],3));
			} else {
				C_temp(i*4+0) = 0;
				C_temp(i*4+1) = 0;
				C_temp(i*4+2) = 0;
				C_temp(i*4+3) = 0;
        // Don't need S_temp
			}
      F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
      F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
      F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
      F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
    }
    
    igl::writeOFF("viewer.off", V_temp, F_temp);
    // Set the viewer data.
    viewer.data.clear();
    viewer.data.set_mesh(V_temp,F_temp);
    viewer.data.set_face_based(true);
		if (_TC.rows() > 0) {
			viewer.data.set_colors(C_temp);
      viewer.core.lighting_factor = 0;
		}
  }


  return false;
}

} // namespace

void TetMeshViewer::viewTetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF) {
  _TV = TV;
  _TT = TT;
  _TF = TF;
  igl::barycenter(_TV,_TT,_B);
  
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down_depth;
  key_down_depth(viewer,'5',0);
  viewer.launch();
}

void TetMeshViewer::viewTetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF,
            const Eigen::VectorXd &TC, bool normalize_cols) {
  _TV = TV;
  _TT = TT;
  _TF = TF;
  _TC = TC;
  igl::barycenter(_TV,_TT,_B);
  
  igl::jet(TC, normalize_cols, _C);
  
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down_depth;
  key_down_depth(viewer,'5',0);
  viewer.launch();
}

void TetMeshViewer::extractOffset(float offset, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

  //igl::marching_cubes(*_C, *_TV, 
}
