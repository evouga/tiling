#include "viewTetMesh.h"

#include <set>

#include <igl/barycenter.h>
#include <igl/bfs_orient.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/jet.h>
#include <igl/viewer/Viewer.h>

#include "offsetSurface.h"
#include "glob_defs.h"

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
  Eigen::VectorXi _TM; // tet markers (contour ids).
  // For the offset surfaces
  Eigen::MatrixXd _OFF_V;
  Eigen::MatrixXi _OFF_F;

  void setValueViewer(double val, igl::viewer::Viewer& viewer) {

    using namespace Eigen;

    std::vector<int> s;

    for (unsigned i=0; i<_TT.rows(); ++i) {
      if (_TC(_TT(i, 0)) <= val && _TC(_TT(i, 1)) <= val &&
          _TC(_TT(i, 2)) <= val && _TC(_TT(i, 3)) <= val) {
        s.push_back(i);
			}
		}

    MatrixXd V_temp(s.size()*4,3);
    MatrixXi F_temp(s.size()*4,3);
		MatrixXd C_temp(s.size()*4,3);

    for (unsigned i=0; i<s.size();++i)
    {
      V_temp.row(i*4+0) = _TV.row(_TT(s[i],0));
      V_temp.row(i*4+1) = _TV.row(_TT(s[i],1));
      V_temp.row(i*4+2) = _TV.row(_TT(s[i],2));
      V_temp.row(i*4+3) = _TV.row(_TT(s[i],3));
      C_temp.row(i*4+0) = _C.row(_TT(s[i],0));
      C_temp.row(i*4+1) = _C.row(_TT(s[i],1));
      C_temp.row(i*4+2) = _C.row(_TT(s[i],2));
      C_temp.row(i*4+3) = _C.row(_TT(s[i],3));
      F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
      F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
      F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
      F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
    }

    // Set the viewer data.
    viewer.data.clear();
    viewer.data.set_mesh(V_temp,F_temp);
    viewer.data.set_face_based(true);
    viewer.data.set_colors(C_temp);
    viewer.core.lighting_factor = 0;
  }

  // This will set the vertices and faces in the viewer,
  // based on the depth. Uses barycenters (stored in _B) to determine
  // what is displayed.
  void setDepthViewer(double depth, igl::viewer::Viewer& viewer) {
    using namespace Eigen;

    Eigen::VectorXd v = _B.col(2).array() - _B.col(2).minCoeff();
    v /= v.col(0).maxCoeff();

    std::vector<int> s;

    for (unsigned i=0; i<v.size(); ++i) {
      if (v(i) <= depth) {
        s.push_back(i);
			}
		}

    MatrixXd V_temp(s.size()*4,3);
    MatrixXi F_temp(s.size()*4,3);
		MatrixXd C_temp(s.size()*4,3);

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
			} else {
				C_temp(i*4+0) = 0;
				C_temp(i*4+1) = 0;
				C_temp(i*4+2) = 0;
				C_temp(i*4+3) = 0;
			}
      F_temp.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
      F_temp.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
      F_temp.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
      F_temp.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
    }

    // Set the viewer data.
    viewer.data.clear();
    viewer.data.set_mesh(V_temp,F_temp);
    viewer.data.set_face_based(true);
		if (_TC.rows() > 0) {
			viewer.data.set_colors(C_temp);
      viewer.core.lighting_factor = 0;
		}
    std::set<int> used;
    for (int i = 0; i < _TM.rows(); i++) {
      if (used.find(_TM(i)) == used.end() && _TM(i) != GLOBAL::nonoriginal_marker) {
        used.insert(_TM(i));
        viewer.data.add_label(_TV.row(i), std::to_string(_TM(i)));
      }
    }
  }

// This function is called every time a keyboard button is pressed
bool key_down_depth(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
  using namespace std;
  using namespace Eigen;

  if (key >= '1' && key <= '9')
  {
    double t = double((key - '1')+1) / 9.0;
    setDepthViewer(t, viewer);
  }


  return false;
}

} // namespace

void TetMeshViewer::viewTetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF) {
  _TV = TV;
  _TT = TT;
  _TF = TF;
  _TC.resize(0);  // Must do or it will suppose there are colors
  igl::barycenter(_TV,_TT,_B);
  
  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down_depth;
  key_down_depth(viewer,'5',0);
  viewer.launch();
}

void TetMeshViewer::viewTetMesh(const Eigen::MatrixXd &TV,
                                const Eigen::MatrixXi &TT,
                                const Eigen::MatrixXi &TF,
                                const Eigen::VectorXd &TC,
                                bool normalize_cols) {
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

void TetMeshViewer::viewTetMesh(const Eigen::MatrixXd &TV,
                                const Eigen::MatrixXi &TT,
                                const Eigen::MatrixXi &TF,
                                const Eigen::VectorXi &TM,
                                const Eigen::VectorXd &TC) {
  _TV = TV;
  _TT = TT;
  _TF = TF;
  _TC = TC;
  _TM = TM;
  igl::barycenter(_TV,_TT,_B);
  igl::jet(TC, true, _C);

  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down_depth;
  key_down_depth(viewer,'5',0);
  viewer.launch();
}

void TetMeshViewer::viewOffsetSurface(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXd &Z) {
  _TV = TV;
  _TT = TT;
  _TF = TF;
  _TC = Z;
  igl::barycenter(_TV,_TT,_B);
  igl::jet(Z, true, _C);
  
  double offset = 0.8;
  OffsetSurface::Triangulation T;
  printf("Generating offset surface...\n");
  OffsetSurface::generateOffsetSurface(TV, TT, Z, offset, _OFF_V, _OFF_F, T);
  igl::viewer::Viewer viewer;
  //viewer.callback_init = [&T, &offset](igl::viewer::Viewer& viewer) {
  printf("Setting up viewer...\n");
  viewer.callback_init = [&](igl::viewer::Viewer& viewer) {
    viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Offset Options");
    viewer.ngui->addVariable("double", offset);
    viewer.ngui->addButton("Redraw",[&offset, &viewer, &T, &TV](){
        //setDepthViewer(offset, viewer); 
        printf("Attempting to draw offset surface at %lf\n", offset);
        OffsetSurface::generateOffsetSurface(T, offset, _OFF_V, _OFF_F);
        viewer.data.clear();
        viewer.data.set_mesh(_OFF_V, _OFF_F);
        });
    viewer.screen->performLayout();
    return false;
  };
  printf("Setting up callback keys "
         "(press 'T' to toggle normal mesh, press 'O' to toggle offset mesh)\n");
  viewer.callback_key_down = [&](igl::viewer::Viewer& viewer, unsigned char key, int modifier) {
      if (key == 'T') {
        setValueViewer(offset, viewer);
      }
      if (key == 'O') {
        viewer.data.clear();
        viewer.data.set_mesh(_OFF_V, _OFF_F);
      }
      return false;
    };
  printf("Success! Now generating mesh\n");
  //OffsetSurface::generateOffsetSurface(TV, Z, 0.5, Voff, Foff);
  //viewer.data.clear();
  //viewer.data.set_mesh(Voff, Foff);
  //viewer.callback_key_down = &key_down_depth;
  //key_down_depth(viewer,'5',0);
  viewer.data.set_mesh(_OFF_V, _OFF_F);
  viewer.launch();
}

void extractOffset(float offset, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

  //igl::marching_cubes(*_C, *_TV, 
}

