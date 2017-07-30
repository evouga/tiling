#ifndef CURVATURE_FLOW_H
#define CURVATURE_FLOW_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

void computeCurvatureFlow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::VectorXi &M,
                          double timestep,
                          Eigen::MatrixXd &Vc);

// Adds extra vertices to the top and bottom to prevent collapsing, then calls
// biharmonic (main function, below) with all the correct inputs.
double biharmonic_new(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                      const Eigen::VectorXi &O,
                      Eigen::MatrixXd &Vc, Eigen::MatrixXi &Fc,
                      Eigen::VectorXi &Mc, std::vector<int> *new_vertices=NULL,
                      bool do_min = false);

// Given some input mesh V,F, produce an output Vc
// Also includes a flag to remove any interior "original" markers, so the
// normals are not messed up.
double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &orig,
                  Eigen::MatrixXd &Vc, bool remove_interior = true);
/**
 * Solve biharmonic for a given mesh specified by V,F. orig contains per-vertex
 * markers, where everything that is fixed is something other than
 * GLOBAL::nonoriginal_marker. Also includes a parameter that has a non-zero value
 * for indices that should have a fixed z-component.
 *
 * @param V,F vertices and faces of mesh
 * @param orig vector of length #V that has a GLOBAL::nonoriginal_marker value for
 *        non-fixed vertices
 * @param fixed_xy a vector of length #V that has a GLOBAL::nonoriginal_marker value
 *        for points whose xy-coordinate should be fixed. This can be the same
 *        as orig.
 * @param L the Laplacian to use.
 * @param out Vc the updated vertices
 * @return the biharmonic energy of the system.
 */
double biharmonic(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::VectorXi &orig, const Eigen::VectorXi &xy_nonfixed,
                  const Eigen::SparseMatrix<double> &L,
                  Eigen::MatrixXd &Vc);

// Same function as above, for viewing.
void biharmonic_view(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::VectorXi &orig,
                Eigen::MatrixXd &Vc,
                bool remove_interior = true);

Eigen::VectorXd biharmonic_energy_per_vertex(const Eigen::MatrixXd &V,
                                             const Eigen::MatrixXi &F,
                                             const std::vector<int> &to_ignore);

// Gives the total biharmonic energy of the system, which includes geodesic
// curvature. Computed as:
//
//   \int a^2 + b^2 dA.
// 
// or
//   4.0 * mean_squared - 2.0 * geodesic;
//
double biharmonic_energy(const Eigen::MatrixXd &V,
                         const Eigen::MatrixXi &F,
                         const std::vector<int> *to_ignore=NULL);

double geodesic_curvature(const Eigen::MatrixXd &V,
                          const Eigen::MatrixXi &F,
                          const std::vector<int> *border_vertices);

#endif // CURVATURE_FLOW_H
