#ifndef TILING_UTILS_H
#define TILING_UTILS_H

#include <vector>
#include <set>

#include <Eigen/Core>

namespace TilingUtils {

struct ConnectedComponent {
  Eigen::MatrixXd V;  // vertices
  Eigen::MatrixXi F;  // faces
  Eigen::VectorXi M;  // markers

  double offsetVal;
  std::set<int> used;

  ConnectedComponent(double offset,
                     const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                     const Eigen::VectorXi &M,
                     const std::set<int> &vertices_used);
};

/**
 * Given a triangle-mesh, will extract all connected components.
 *
 * Arguments:
 *  @param vertex_index - the vertex to dfs from.
 *  @param visited - the vertices visited.
 *  @param vertices_used - the vertices in the component.
 *  @param F - faces.
 */
void dfs(int vertex_index, std::set<int>& visited,
         std::set<int>& vertices_used,
         const Eigen::MatrixXi &F);

/**
 * Given a triangle-mesh, will extract all connected components.
 *
 * Arguments:
 *  @param V - vertices.
 *  @param F - faces.
 *  @param O - markers.
 *  @param offset - offset value used to generate the mesh.
 *
 * Returns:
 *  A vector of connected components.
 */
std::vector<ConnectedComponent> getConnectedComponents(const Eigen::MatrixXd &V,
                                                       const Eigen::MatrixXi &F,
                                                       const Eigen::VectorXi &O,
                                                       double offset);

/**
 * Does a depth-first search on the entire mesh to gather contour points.
 *
 * Arguments:
 *  @param vertex_index - the vertex to dfs from.
 *  @param visited - the vertices visited.
 *  @param vertices_used - the vertices in the component.
 *  @param F - faces.
 *  @param O - original markers.
 */
void addVerticesToContour(int vertex_index,
                          std::set<int>& visited, std::set<int> &vertices_used,
                          const Eigen::MatrixXi &F, const Eigen::VectorXi &O);

/**
 * Given a tile, returns the original contours.
 *
 * Arguments:
 *  @param F - faces.
 *  @param O - markers.
 *
 * Returns:
 *  A vector of sets that contain the vertices in the contour.
 */
std::vector<std::set<int> > getContourVertices(const Eigen::MatrixXi &F,
                                               const Eigen::VectorXi &O);

/**
 * Returns the contours used from a tile.
 *
 * Arguments:
 *  @param component - tile to extract original contours.
 *  @param contours - original contours.
 *  @param V - vertices corresponding to the contours.
 *
 * Returns:
 *  Set of contours used.
 */
std::set<int> getContoursUsed(const ConnectedComponent &component,
                              std::vector<std::set<int> > &contours,
                              const Eigen::MatrixXd &V);

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
