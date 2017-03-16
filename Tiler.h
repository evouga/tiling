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
  Tile *bottom_parent;

  bool is_last = false;

  // Used for get upper connected, since it is recursive.
  bool cached_upper_connected = false;
  std::vector<std::set<int> > upper_connected;

  // Contours to tile.
  const std::set<int> upper;

  // Configuration of connections between contours and parent.
  std::vector<Component*> components;

  // Constructor for the starting level (no parents).
  Tile(const std::set<int> &upper_contours);

  // Constructor for tiles that have parent constraints.
  Tile(const std::set<int> &upper_contours,
       const std::vector<std::set<int> > &components,
       Tile *parent, bool last);

  bool isTopTile() { return this->is_last; }
  bool isValid();

  // Forms a graph depicting the current tile (includes parent connectivity.
  std::map<int, std::vector<int> > getGraph();

  // Each set says if the contours were ever connected (including parents).
  std::vector<std::set<int> > getUpperConnected();

  // Gets a full tile to the very bottom most parent.
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
 *  @param is_last - additional topology constraints on the last tile.
 */
std::vector<Tile*> generateTiles(const std::set<int> &upper, Tile *parent,
                                 bool is_last);

/**
 * Populates all valid configurations of connections.
 * Use this one if you want to pass in the sets  used in the set cover.
 *
 * Arguments:
 *  @param upper - upper contours to be covered.
 *  @param parent - can be NULL, the tile below.
 *  @param components - pass in the possible sets in the cover.
 *  @param is_last - additional topology constraints on the last tile.
 */
std::vector<Tile*> generateTiles(const std::set<int> &upper, Tile *parent,
                                 const std::vector<std::set<int> > &components,
                                 bool is_last);

/**
 * Returns various levels' components of a tile.
 *
 * Arguments:
 *  @param tile - tile to get mesh on.
 *  @param num_tiles - 0 for current layer, 1 for two layers, -1 for all layers.
 */
std::vector<Component*> getTileComponents (Tile *tile, int num_tiles);

} // namespace Tiler

#endif
