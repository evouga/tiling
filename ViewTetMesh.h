#include <igl/viewer/Viewer.h>
#include <igl/barycenter.h>

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);
void loadTetMesh (Eigen::MatrixXd TV1, Eigen::MatrixXi TT1, Eigen::MatrixXi TF1);
void loadTetMesh (const Eigen::MatrixXd &TV1, const Eigen::MatrixXi &TT1, const Eigen::MatrixXi &TF1,
									const Eigen::MatrixXd &C);
