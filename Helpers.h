#ifndef HELPERS_H
#define HELPERS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Helpers {

bool is_edge_manifold(const Eigen::MatrixXi &F, Eigen::VectorXi &M);

// For making the mesh manifold.
void extractManifoldPatch(
    Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O,
    int minFaces=5, // must have more than one tet (4 faces)
    bool changeMarkers=true); 

void removeDuplicates(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O);

Eigen::VectorXi getFaceMarkers(const Eigen::MatrixXi& F, const Eigen::VectorXi& M);

// TODO(bradyz): Create a tri-mesh data structure.
void tetrahedralize(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                    const Eigen::VectorXi &VM, const Eigen::VectorXi &FM,
                    Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
                    Eigen::MatrixXi &TF, Eigen::VectorXi &TO);

void triangulate(const Eigen::MatrixXd &P, const Eigen::MatrixXi &E,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void viewTriMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                 const Eigen::VectorXi &O);

void offsetViewer(const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TT,
                  const Eigen::VectorXd &VH);

void combineMesh(const std::vector<Eigen::MatrixXd> &Vs,
                 const std::vector<Eigen::MatrixXi> &Fs,
                 const std::vector<Eigen::VectorXi> &Ms,
                 Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &M);

void extractShell(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1,
                  const Eigen::VectorXi &M1,
                  Eigen::MatrixXd &V2, Eigen::MatrixXi &F2, Eigen::VectorXi &M2);

void collapseSmallTriangles(const Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                            double eps=1e-8);

void removeUnreferenced(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O);

bool isMeshOkay(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                double eps=1e-8);

bool sparseMatrixHasNaN(const Eigen::SparseMatrix<double> &A);

} // namespace Helpers

#endif // HELPERS_H
