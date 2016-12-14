#ifndef TILING_UTILS_H
#define TILING_UTILS_H

#include <vector>

#include <Eigen/Core>

namespace TilingUtils {
struct ConnectedComponent {
  Eigen::MatrixXd V;  // vertices
  Eigen::MatrixXd F;  // faces
  Eigen::VectorXi O;  // original?

  double offsetVal;
};

/**
 * Will extract all possible Connected Components based off of heat values of
 * a given tet-mesh (3d).
 *
 * Arguments are:
 * @param TV - vertices of tet-mesh
 * @param TF - faces of tet-mesh
 * @param TT - tets of tet-mesh
 * @param TO - original markers list
 * @param H - heat values for each vertex in TV
 */
std::vector<ConnectedComponent> allPossibleTiles(
    const Eigen::MatrixXd &TV, const Eigen::MatrixXi &TF, const Eigen::MatrixXi &TT,
    const Eigen::VectorXi &TO, const Eigen::VectorXd &H);

} // namespace TilingUtils
#endif //TILING_UTILS_H
