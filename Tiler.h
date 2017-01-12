#ifndef TILER_H
#define TILER_H

#include <set>
#include <vector>

#include <Eigen/Core>

namespace Tiler {

// Consists of one to many contours.
struct Component {
  const std::set<int> contours_used;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::VectorXi M;

  Component(const std::set<int> contours) : contours_used(contours) { }
};

// A tile is a set of components that cover all upper and lower contours.
struct Tile {
  // The tile that comes right before the current one.
  const Tile *bottom_parent;

  // TODO(bradyz): Needs a fixup.
  bool isLast = false;

  // Contours to tile.
  const std::set<int> upper;

  // Configuration of connections between contours and parent.
  std::vector<Component*> components;

  // First tile constructor.
  Tile(const std::set<int> &upper_contours);

  // Tile with parent constructor.
  Tile(const std::set<int> &upper_contours,
       const std::vector<std::set<int> > &components,
       const Tile *parent);

  bool isTopTile() { return this->isLast; }
  bool isValid() const;

  // TODO(bradyz): get rid of this function.
  std::vector<std::set<int> > getUpperConnected() const;

  std::vector<Component*> getAllComponents() const;

  ~Tile();
};

// Ostream overloads for debugging.
std::ostream& operator<< (std::ostream &os, const Component &comp);
std::ostream& operator<< (std::ostream &os, const Tile &tile);

/**
 * Populates all valid configurations of connections.
 *
 * Arguments:
 *  @param upper - upper contours to be covered.
 *  @param parent - can be NULL, the tile below.
 *  @param isLast - additional topology constraints on the last tile.
 */
std::vector<Tile*> generateTiles(const std::set<int> &upper, const Tile *parent,
                                 bool isLast);

/**
 * Populates all valid configurations of connections.
 * Use this one if you want to pass in the sets  used in the set cover.
 *
 * Arguments:
 *  @param upper - upper contours to be covered.
 *  @param parent - can be NULL, the tile below.
 *  @param components - pass in the possible sets in the cover.
 *  @param isLast - additional topology constraints on the last tile.
 */
std::vector<Tile*> generateTiles(const std::set<int> &upper, const Tile *parent,
                                 const std::vector<std::set<int> > &components,
                                 bool isLast);

/**
 * View a full tile.
 *
 * Arguments:
 *  @param tile - tile to be viewed.
 */
void viewTile(const Tile *tile);

/**
 * Combine all components and extract the surface of a tile.
 *
 * Arguments:
 *  @param tile - tile to get mesh on.
 *  @param V - output vertices.
 *  @param F - output faces.
 *  @param O - output markers.
 */
void getTileMesh (Tile *tile,
                  Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &O);

} // namespace Tiler

#endif
