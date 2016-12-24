// This file reads from standard in and creates a stack of tiles.
// Input format is -
// N total tiles.
// for i from 1 to N
//   contours covered
//   M number of connected components.
//   for j from 1 to M
//     contours used in connected component i.
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>

#include "../Tiler.h"

using namespace std;
using namespace Tiler;

int levels;

Tile *readTile(Tile *parent) {
  int n;
  string buffer;

  getline(cin, buffer);

  set<int> upper;
  stringstream ss(buffer);
  while (ss >> n)
    upper.insert(n);

  cin >> n;
  getline(cin, buffer);

  vector<Component*> components;
  for (int i = 0; i < n; i++) {
    getline(cin, buffer);
    ss = stringstream(buffer);

    set<int> comp;
    int contour_id;
    while (ss >> contour_id)
      comp.insert(contour_id);
    components.push_back(new Component(comp));
  }

  Tile *tile = new Tile(upper);
  tile->components = components;

  return tile;
}

int main() {
  cin >> levels;
  string buffer;
  getline(cin, buffer);

  // Parse in all the tiles.
  Tile **tiles = new Tile*[levels];
  for (int i = levels-1; i >= 0; i--) {
    tiles[i] = readTile(NULL);
    if (i != levels-1)
      tiles[i+1]->bottom_parent = tiles[i];
  }

  cout << *tiles[levels-1] << endl;
  cout << "is valid: " << tiles[levels-1]->isValid() << endl;
}
