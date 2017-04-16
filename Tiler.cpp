#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Tiler.h"
#include "Helpers.h"

using namespace std;

// Helper functions, bit shifting stuff and exact set cover.
namespace {

// Having a bunch of unsigned ints make lines long.
typedef unsigned int uint;

bool isPowerOfTwo(uint x) {
  return (x != 0 && !(x & (x - 1)));
}

bool getNthBit(uint x, uint n) {
  return ((x >> n) & 1);
}

uint setNthBit(uint x, uint n) {
  return (x | (1 << n));
}

string binaryRepresentation(uint x, uint binary_width) {
  char buffer[binary_width];
  for (int i = binary_width-1; i >= 0; i--)
    buffer[i] = (getNthBit(x, i) ? '1' : '0');
  return string(buffer);
}

set<int> setRepresentation(uint x, uint binary_width, const vector<int> &data) {
  set<int> result;
  for (int i = binary_width-1; i >= 0; i--)
    if (getNthBit(x, i))
      result.insert(data[i]);
  return result;
}

template <typename T>
string matrixToString(T **matrix, size_t rows, size_t cols) {
  stringstream ss;
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      ss << matrix[i][j] << " ";
    }
    if (i != rows-1)
      ss << endl;
  }
  return ss.str();
}

string graphToString(const map<int, vector<int> > &graph) {
  stringstream ss;
  for (const auto &it : graph) {
    int u = it.first;
    ss << u << ": ";
    for (int v : it.second)
      ss << v << " ";
    ss << endl;
  }
  return ss.str();
}

set<int> set_union(const set<int> &S, const set<int> &T) {
  set<int> both;
  for (int x : S)
    both.insert(x);
  for (int x : T)
    both.insert(x);
  return both;
}

void exactSetCoverUtility(bool **A, size_t n, size_t m,
                          set<uint> &deleted_rows, set<uint> &deleted_cols,
                          set<uint> &used, set<set<uint> > &set_covers) {
  // Found a valid cover.
  if (deleted_cols.size() == m) {
    set_covers.insert(used);
    return;
  }

  for (uint i = 0; i < n; i++) {
    if (deleted_rows.find(i) != deleted_rows.end())
      continue;
    bool has_extra = false;
    for (uint j : deleted_cols) {
      has_extra |= A[i][j];
    }
    if (has_extra)
      continue;

    // Update current result.
    used.insert(i);
    deleted_rows.insert(i);
    for (uint j = 0; j < m; j++) {
      if (A[i][j])
        deleted_cols.insert(j);
    }
    // Recurse.
    exactSetCoverUtility(A, n, m,
                         deleted_rows, deleted_cols, used, set_covers);
    // Backtrack.
    used.erase(i);
    deleted_rows.erase(i);
    for (uint j = 0; j < m; j++) {
      if (A[i][j])
        deleted_cols.erase(j);
    }
  }
}

// Algorithm X, Dancing links (really cool names, too bad it's NP-hard).
set<set<uint> > exactSetCover(const vector<uint> &subsets, uint m) {
  uint n = subsets.size();

  // A(i, j) is true if subset i covers element j.
  bool **matrix = new bool*[n];
  for (size_t i = 0; i < n; i++) {
    matrix[i] = new bool[m];
    for (size_t j = 0; j < m; j++)
      matrix[i][j] = getNthBit(subsets[i], j);
  }

  // Populated by utility function.
  set<set<uint> > set_covers;
  set<uint> rows_used;
  set<uint> cols_used;
  set<uint> used;
  exactSetCoverUtility(matrix, n, m, rows_used, cols_used, used, set_covers);

  // Cleanup.
  for (size_t i = 0; i < n; i++)
    delete matrix[i];
  delete[] matrix;

  return set_covers;
}

template <typename T>
vector<vector<set<T> > > exactSetCover(set<T> to_cover, vector<set<T> > subsets) {
  // Generate element ids {0, 1, ..., n}.
  uint current_size = 0;
  map<T, uint> to_binary;
  for (set<T> &subset : subsets)
    for (T element : subset)
      if (to_binary.find(element) == to_binary.end())
        to_binary[element] = current_size++;
  for (T element : to_cover)
    if (to_binary.find(element) == to_binary.end())
      to_binary[element] = current_size++;

  // Convert the subsets to their binary representation.
  // If we have 3 total elements:
  // 01 => 011
  // 2  => 100
  vector<uint> subsets_binary;
  for (set<T> &subset : subsets) {
    uint subset_binary = 0;
    for (T element : subset)
      subset_binary = setNthBit(subset_binary, to_binary[element]);
    subsets_binary.push_back(subset_binary);
  }

  // Algorithm X.
  set<set<uint> > covers = exactSetCover(subsets_binary, to_binary.size());

  // Convert the binary covers into their T representation.
  vector<vector<set<T> > > result;
  for (set<uint> solution : covers) {
    vector<set<T> > cover;
    for (uint i : solution)
      cover.push_back(subsets[i]);
    result.push_back(cover);
  }

  return result;
}

bool isTree(int start, const map<int, vector<int> > &graph, set<int> &visited) {
  visited.insert(start);

  vector<int> the_stack;
  the_stack.push_back(start);

  int nodes_visited = 0;
  int edges_visited = 0;

  // DFS to find the connected component.
  while (the_stack.size() > 0) {
    int u = the_stack.back();
    the_stack.pop_back();

    // Note: double counting edges, need to half later.
    nodes_visited++;
    edges_visited += graph.at(u).size();

    // Visit neighbors.
    for (int v : graph.at(u)) {
      if (visited.find(v) != visited.end())
        continue;
      visited.insert(v);
      the_stack.push_back(v);
    }
  }

  edges_visited /= 2;

  return (edges_visited == (nodes_visited - 1));
}

bool isForest(const map<int, vector<int> > &graph) {
  set<int> visited;
  for (const auto &it : graph) {
    int u = it.first;
    if (visited.find(u) != visited.end())
     continue;
    else if (!isTree(u, graph, visited))
      return false;
  }
  return true;
}

bool isConnected(const set<int> &upper, const set<int> &lower,
                 const vector<vector<int> > &connectedComponents) {
  for (int u : lower) {
    for (const vector<int> &component : connectedComponents) {
      if (find(component.begin(), component.end(), u) != component.end()) {
        bool valid = false;
        for (int v : component)
          valid |= (upper.find(v) != upper.end());
        if (!valid)
          return false;
        break;
      }
    }
  }
  return true;
}

void dfs(int u, set<int> &visited, const map<int, vector<int> > &graph,
         vector<int> & component) {
  visited.insert(u);
  component.push_back(u);
  for (int v : graph.at(u)) {
    if (visited.find(v) == visited.end())
      dfs(v, visited, graph, component);
  }
}

vector<vector<int> > getConnectedComponents(const map<int, vector<int> > &graph) {
  set<int> visited;
  vector<vector<int> > result;
  for (const auto &it : graph) {
    int u = it.first;
    if (visited.find(u) == visited.end()) {
      vector<int> component;
      dfs(u, visited, graph, component);
      if (component.size() > 0)
        result.push_back(component);
    }
  }
  return result;
}

// The following optimizations can be done if the inputs are (0, 12, 3).
// - 0 cannot exist alone, reject 0001.
// - 3 cannot exist alone, reject 1000.
// - 12 cannot exist in the same subset, reject all of the form *11*.
//
// Current runtime is exponential.
vector<set<int> > generatePossibleSubsets(const set<int> &entire_set,
                                          const vector<set<int> > &previous) {
  vector<int> entire_set_vector(entire_set.begin(), entire_set.end());
  int total = entire_set.size();
  vector<set<int> > possible;
  for (unsigned int counter = 1; counter < (2 << total); counter++) {
    // Optimization 1: Ignore single elements.
    if (isPowerOfTwo(counter))
      continue;

    // Optimization 2: Ignore loops.
    bool can_use = true;
    for (int i = 0; i < total && can_use; i++) {
      if (!getNthBit(counter, i))
        continue;
      for (int j = i+1; j < total && can_use; j++) {
        if (!getNthBit(counter, j))
          continue;
        for (const set<int> &connected : previous) {
          bool used1 = (connected.find(i) != connected.end());
          bool used2 = (connected.find(j) != connected.end());
          if (used1 && used2)
            can_use = false;
        }
      }
    }

    if (can_use)
      possible.push_back(setRepresentation(counter, total, entire_set_vector));
  }
  return possible;
}

} // anonymous namespace

namespace Tiler {

Tile::Tile(const set<int> &upper_contours) :
    bottom_parent(NULL),
    is_last(false),
    upper(upper_contours) {
  for (int x : upper_contours) {
    set<int> component;
    component.insert(x);
    this->components.push_back(new Component(component));
  }
}

Tile::Tile(const set<int> &upper_contours, const vector<set<int> > &components,
           Tile *parent, bool last) :
    bottom_parent(parent),
    is_last(last),
    upper(upper_contours) {
  for (const set<int> &contours : components)
    this->components.push_back(new Component(contours));
}

// TODO(bradyz): The bug in cycles is here.
vector<set<int> > Tile::getUpperConnected() {
  // The first slice is all disjoint.
  if (this->bottom_parent == NULL) {
    for (int x : this->upper) {
      set<int> component;
      component.insert(x);
      this->upper_connected.push_back(component);
    }

    // Save for future use (I think the term is "thunking").
    this->cached_upper_connected = true;
    return this->upper_connected;
  }

  // Save some recursive calls.
  if (this->cached_upper_connected)
    return this->upper_connected;

  map<int, vector<int> > graph = this->getGraph();
  vector<vector<int> > connectedComponents = getConnectedComponents(graph);

  for (const vector<int> &component : connectedComponents) {
    set<int> component_upper;

    for (int u : component) {
      if (this->upper.find(u) != this->upper.end())
        component_upper.insert(u);
    }

    if (component_upper.size() > 0)
      this->upper_connected.push_back(component_upper);
  }

  // Save the result.
  this->cached_upper_connected = true;
  return this->upper_connected;
}

vector<Component*> Tile::getAllComponents() const {
  vector<Component*> result;

  if (this->bottom_parent != NULL) {
    result = this->bottom_parent->getAllComponents();

    for (Component *component : this->components)
      result.push_back(component);
  }

  return result;
}

ostream& operator<< (ostream &os, const Component &comp) {
  for (int x : comp.contours_used)
    os << x << " ";
  return os;
}

ostream& operator<< (ostream &os, const Tile &tile) {
  os << tile.components.size() << " connected components covering - ";
  for (int x : tile.upper)
    os << x << " ";
  os << endl;
  for (int i = 0; i < tile.components.size(); i++) {
    os << *tile.components[i];
    if (i != (tile.components.size() - 1))
      os << endl;
  }
  if (tile.bottom_parent != NULL)
    return os << endl << *tile.bottom_parent;
  return os;
}

bool Tile::isValid() {
  if (this->bottom_parent == NULL)
    return true;

  map<int, vector<int> > graph = this->getGraph();
  vector<vector<int> > connectedComponents = getConnectedComponents(graph);

  bool is_valid = true;

  // The contours must not form any cycles.
  is_valid &= isForest(graph);

  // The bottom and top must be connected in at least one point.
  is_valid &= isConnected(this->upper, bottom_parent->upper, connectedComponents);

  // When the tile is done, there must be only one connected component.
  is_valid &= (!this->is_last || (connectedComponents.size() == 1));

  // cout << "graph: " << endl;
  // cout << graphToString(graph) << endl;

  return is_valid;
}

map<int, vector<int> > Tile::getGraph() {
  // Outputs an undirected graph that points from bottom to top.
  // For example, if the bottom has {1}, and the top has {2, 3},
  // and there is one connected component {1, 2, 3},
  // the graph would look like the following -
  //
  // 2   3
  //  \ /
  // 1000
  //   |
  //   1

  map<int, vector<int> > graph;

  // There are no edges in the bottom most tile.
  if (this->bottom_parent == NULL) {
    for (int u : this->upper)
      graph[u] = vector<int>();
    return graph;
  }

  Tile* parent = this->bottom_parent;

  // Add new edges from components in tile.
  // 1*** means connected in the current layer.
  int dummy_node = 1000;
  for (const Component *component : this->components) {
    if (component->contours_used.size() == 1) {
      int x = *component->contours_used.begin();
      if (parent->upper.find(x) != parent->upper.end()) {
        graph[x] = vector<int>();
        continue;
      }
    }

    for (int u : component->contours_used) {
      graph[u].push_back(dummy_node);
      graph[dummy_node].push_back(u);
    }
    dummy_node++;
  }

  // 2*** means at one point connected.
  dummy_node = 2000;
  for (const set<int> &connected : this->bottom_parent->getUpperConnected()) {
    if (connected.size() == 1)
      continue;
    for (int u : connected) {
      graph[u].push_back(dummy_node);
      graph[dummy_node].push_back(u);
    }
    dummy_node++;
  }

  return graph;
}

Tile::~Tile() {
  for (Component* comp : this->components)
    delete comp;
}

vector<Tile*> generateTiles(const set<int> &upper, Tile *parent,
                            bool is_last) {
  vector<Tile*> result;

  // Bottom tile.
  if (parent == NULL) {
    result.push_back(new Tile(upper));
    return result;
  }

  // Contours to cover.
  set<int> all_contours = set_union(upper, parent->upper);

  // These are used to generate the set cover.
  vector<set<int> > previous = parent->getUpperConnected();
  vector<set<int> > components = generatePossibleSubsets(all_contours, previous);

  // Generate possible tilings.
  vector<vector<set<int> > > possible = exactSetCover(all_contours, components);

  // Go through all possible tilings and find ones with no loops.
  for (const vector<set<int> > &tiling : possible) {
    Tile *tile = new Tile(upper, tiling, parent, is_last);

    // Get rid of the tile if it's no good.
    if (tile->isValid())
      result.push_back(tile);
    else
      delete tile;
  }

  return result;
}

// Use this function to pass in valid components used in the cover.
vector<Tile*> generateTiles(const set<int> &upper, Tile *parent,
                            const vector<set<int> > &components, bool is_last) {
  vector<Tile*> result;

  // Bottom tile.
  if (parent == NULL) {
    result.push_back(new Tile(upper));
    return result;
  }

  // Contours to cover.
  set<int> all_contours = set_union(upper, parent->upper);

  // Generate possible tilings.
  vector<vector<set<int> > > possible = exactSetCover(all_contours, components);

  // Go through all possible tilings and find ones with no loops.
  for (const vector<set<int> > &tiling : possible) {
    Tile *tile = new Tile(upper, tiling, parent, is_last);
    tile->is_last = is_last;

    // Get rid of the tile if it's no good.
    if (tile->isValid())
      result.push_back(tile);
    else
      delete tile;
  }

  return result;
}

// Extracts the surface of the current level of the given tile.
vector<Component*> getTileComponents (Tile *tile,  int num_tiles) {
  vector<Component*> components;

  if (num_tiles == 0) {
    components = tile->getAllComponents();
  }
  else if (num_tiles == 1) {
    components = tile->components;
  }
  else if (num_tiles == 2) {
    // Go through current level.
    for (Component *component : tile->components)
      components.push_back(component);

    // Go through previous level.
    if (tile->bottom_parent != NULL &&
        tile->bottom_parent->bottom_parent != NULL) {
      for (Component *component : tile->bottom_parent->components)
        components.push_back(component);
    }
  }

  return components;
}

} // Namespace Tiler
