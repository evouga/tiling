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
  map<T, uint> to_binary;
  for (set<T> &subset : subsets)
    for (T element : subset)
      if (to_binary.find(element) == to_binary.end())
        to_binary[element] = to_binary.size();
  for (T element : to_cover)
    if (to_binary.find(element) == to_binary.end())
      to_binary[element] = to_binary.size();

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
    if (visited.find(u) == visited.end() && !isTree(u, graph, visited))
      return false;
  }
  return true;
}

bool isConnected(const set<int> &lower, const set<int> &upper,
                 const map<int, vector<int> > &graph) {
  for (const auto &it : graph) {
    const vector<int> &can_reach = it.second;
    for (int j = 0; j < can_reach.size(); j++) {
      for (int k = j+1; k < can_reach.size(); k++) {
        int u = can_reach[j];
        int v = can_reach[k];
        bool u_lower = (lower.find(u) != lower.end());
        bool v_lower = (lower.find(v) != lower.end());
        if (u_lower ^ v_lower)
          return true;
      }
    }
  }
  return false;
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

vector<set<int> > getCorrespondingUpper(const set<int> &upper,
                                        const vector<vector<int> > &components) {
  vector<set<int> > result;
  for (const vector<int> &component : components) {
    set<int> component_upper;
    for (int contour : component) {
      if (upper.find(contour) != upper.end())
        component_upper.insert(contour);
    }
    if (component_upper.size() > 0)
      result.push_back(component_upper);
  }
  return result;
}


} // anonymous namespace

namespace Tiler {

vector<set<int> > generateInitialRequirements(int number_contours) {
  vector<set<int> > result;
  for (int i = 0; i < number_contours; i++) {
    set<int> tmp;
    tmp.insert(i);
    result.push_back(tmp);
  }
  return result;
}

void generateTopAndBottom(int bot_count, int top_count, int offset,
                          set<int> &bot_ids, set<int> &top_ids,
                          set<int> &both) {
  bot_ids.clear();
  top_ids.clear();
  both.clear();

  // Bottom ids go from [0, bot_count).
  int current = offset;
  for (int i = 0; i < bot_count; i++) {
    both.insert(current);
    bot_ids.insert(current);
    current++;
  }

  // Top ids go from [bot_count, bot_count + top_count).
  for (int i = 0; i < top_count; i++) {
    both.insert(current);
    top_ids.insert(current);
    current++;
  }
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

bool isValidTiling(const vector<set<int> > &tile,
                   const set<int> &lower,
                   const set<int> &upper,
                   const vector<set<int> > &previous,
                   bool is_last, vector<vector<int> > &connectedComponents) {
  map<int, vector<int> > graph;

  // Add edges from previous connected components.
  int dummy_node = *upper.rbegin() + 1;
  for (const set<int> &connected : previous) {
    for (int u : connected) {
      graph[u].push_back(dummy_node);
      graph[dummy_node].push_back(u);
    }
    dummy_node++;
  }

  // Add new edges from components in tile.
  for (const set<int> &component : tile) {
    for (int u : component) {
      graph[u].push_back(dummy_node);
      graph[dummy_node].push_back(u);
    }
    dummy_node++;
  }

  connectedComponents = getConnectedComponents(graph);

  bool is_valid = true;
  is_valid &= isForest(graph);
  is_valid &= isConnected(lower, upper, graph);
  is_valid &= (!is_last || (connectedComponents.size() == 1));

  if (is_valid) {
    // cout << "test: " << endl;
    // for (auto x : tile) {
    //   for (auto y : x)
    //     cout << y << " ";
    //   cout << endl;
    // }
    // cout << graphToString(graph) << endl;
  }

  return is_valid;
}

void generateTiles(const set<int> &lower, const set<int> &upper,
                   const set<int> &all_contours,
                   const vector<set<int> > &previous,
                   vector<vector<set<int> > > &result,
                   vector<vector<set<int> > > &upper_used,
                   bool is_last) {
  result.clear();
  upper_used.clear();

  size_t total_nodes = all_contours.size();

  // These are used to generate the set cover.
  vector<set<int> > components = generatePossibleSubsets(all_contours, previous);

  // Generate possible tilings.
  vector<vector<set<int> > > possible = exactSetCover(all_contours, components);

  // Go through all possible tilings and find ones with no loops.
  for (int i = 0; i < possible.size(); i++) {
    const vector<set<int> > &tile = possible[i];
    vector<vector<int> > connected;
    if (isValidTiling(tile, lower, upper, previous, is_last, connected)) {
      result.push_back(tile);
      upper_used.push_back(getCorrespondingUpper(upper, connected));
    }
  }
}

} // namespace Tiler
