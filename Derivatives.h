#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace derivatives {
/*
 * Will decompose L into its Hodge-star definition, which is:
 *
 *   L = D^T \star D
 *
 * where L is the |V|x|V| cotangent Laplacian, D is |V|x|E|, and
 * \star is |E|x|E|. Note that each edge in E is directed (i.e.
 * |E| = |F| * 3, whereas for undirected edges, |E| != |F|*3.
 
 * @param V - vertices of mesh
 * @param F - faces of mesh
 * @param OUT D, star
 */
void decompose_L(const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &V,
                 const Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &F,
                 Eigen::SparseMatrix<double> &D, Eigen::SparseMatrix<double> &star);


/**
 * Use this function for the edges required in these methods.
 */
void buildEdges(const Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &F,
                Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &E);

/**
 * Will compute the geodesic curvature of the given mesh; if dGeodesicCurvature is
 * also supplied, will provide the per-vertex derivative of GeodesicCurvature.
 *
 * If border_vertices are not known, set border_vertices to NULL and they will
 * be computed.
 */
double geodesicCurvature(const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &V,
                         const Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &F,
                         const Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &E,
                         const std::vector<int> *border_vertices,
                         Eigen::Matrix<double, -1, -1, Eigen::RowMajor> *dGeodesicCurvature /* optional */);
/**
 * Helper without edges.
 */
double geodesicCurvature(const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> &V,
                         const Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &F,
                         const std::vector<int> *border_vertices,
                         Eigen::Matrix<double, -1, -1, Eigen::RowMajor> *dGeodesicCurvature /* optional */);

} // namespace derivatives

#endif // DERIVATIVES_H
