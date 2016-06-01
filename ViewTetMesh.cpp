#include <igl/viewer/Viewer.h>
#include <igl/barycenter.h>

// Input polygon
Eigen::MatrixXd _V;
Eigen::MatrixXi _F;
Eigen::MatrixXd _B;

// Tetrahedralized interior
Eigen::MatrixXd _TV;
Eigen::MatrixXi _TT;
Eigen::MatrixXi _TF;
Eigen::MatrixXd _C;

// This function is called every time a keyboard button is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
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
		vector<double> s_col;

    for (unsigned i=0; i<v.size(); ++i) {
      if (v(i) <= t) {
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
			if (_C.rows() > 0) {
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
		/**
		 * This is the ideal method for computing the new faces.
		 * However, tetgen doesn't add any additional faces, so it's
		 * we actually need to do the above method. Which doesn't
		 * work with colors (all the way)
		 *
		 * I'm leaving this section commented out.
		 */
		// Count up the number of faces
		int num_faces = 0;
		for (unsigned int i = 0; i < _TF.rows(); ++i) {
			bool inside = true;
			for (unsigned int j = 0; j < 3; ++j) {
				if (_TV(_TF(i,j),2) > t_abs) {
					inside = false;
					break;
				}
			}
			if (inside) ++num_faces;
		}
    MatrixXi F_temp2(num_faces,3);
		int offset = 0;
		for (unsigned int i = 0; i < _TF.rows(); ++i) {
			bool inside = true;
			for (unsigned int j = 0; j < 3; ++j) {
				if (_TV(_TF(i,j),2) > t_abs) {
					inside = false;
					break;
				}
			}
			if (inside) {
				F_temp2.row(offset++) = _TF.row(i);
			}
		}

		// Count up the number of tets
		int num_tets = 0;
		for (unsigned int i = 0; i < _TT.rows(); ++i) {
			bool inside = true;
			for (unsigned int j = 0; j < 4; ++j) {
				if (_TV(_TT(i,j),2) > t_abs) {
					inside = false;
					break;
				}
			}
			if (inside) ++num_tets;
		}
    MatrixXi T_temp(num_tets,4);
		offset = 0;
		for (int i = 0; i < _TT.rows(); ++i) {
			bool inside = true;
			for (unsigned int j = 0; j < 4; ++j) {
				if (_TV(_TT(i,j),2) > t_abs) {
					inside = false;
					break;
				}
			}
			if (inside) {
				T_temp.row(offset++) = _TT.row(i);
			}
		}
    viewer.data.clear();
    viewer.data.set_mesh(_TV,F_temp2);
		viewer.core.show_lines = false;
		if (_C.rows() > 0) {
			viewer.data.set_colors(_C);
		}
		/*
		 * And following is the previous method for displaying color.
		 */
		/*
    viewer.data.clear();
    viewer.data.set_mesh(V_temp,F_temp);
    viewer.data.set_face_based(true);
		viewer.core.show_lines = false;
		if (_C.rows() > 0) {
			viewer.data.set_colors(C_temp);
		}
		*/
  }


  return false;
}

void loadTetMesh (const Eigen::MatrixXd &TV1, const Eigen::MatrixXi &TT1, const Eigen::MatrixXi &TF1,
									const Eigen::MatrixXd &C) {
  using namespace Eigen;
  using namespace std;

  _TV = TV1;
  _TT = TT1;
  _TF = TF1;
	_C = C;

  igl::barycenter(_TV,_TT,_B);

  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  key_down(viewer,'5',0);
	//viewer.data.set_colors(C);
  viewer.launch();
}

void loadTetMesh (Eigen::MatrixXd TV1, Eigen::MatrixXi TT1, Eigen::MatrixXi TF1) {
  using namespace Eigen;
  using namespace std;

  _TV = TV1;
  _TT = TT1;
  _TF = TF1;

  igl::barycenter(_TV,_TT,_B);

  igl::viewer::Viewer viewer;
  viewer.callback_key_down = &key_down;
  key_down(viewer,'5',0);
  viewer.launch();
}
