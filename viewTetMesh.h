#ifndef VIEW_TET_MESH_H
#define VIEW_TET_MESH_H

#include <Eigen/Core>

namespace TetMeshViewer {

void viewTetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF);
void viewTetMesh(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT, const Eigen::MatrixXi &TF,
            const Eigen::VectorXd &TC, bool normalize_cols);

void viewTetMesh(const Eigen::MatrixXd &TV,
                                const Eigen::MatrixXi &TT,
                                const Eigen::MatrixXi &TF,
                                const Eigen::VectorXi &TM,
                                const Eigen::VectorXd &TC);

// Launch viewer to view offset surfaces
void viewOffsetSurface(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXd &Z);

void extractOffset(float offset, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
} // namespace TetMeshViewer

#endif // VIEW_TET_MESH_H
