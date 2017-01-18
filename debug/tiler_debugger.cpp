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

Tile *readTile(Tile *parent, bool is_last) {
  int n;
  string buffer;

  getline(cin, buffer);
  stringstream ss(buffer);

  set<int> upper;
  while (ss >> n)
    upper.insert(n);

  cin >> n;
  getline(cin, buffer);

  vector<set<int> > components;
  for (int i = 0; i < n; i++) {
    getline(cin, buffer);
    ss = stringstream(buffer);

    set<int> comp;
    int contour_id;
    while (ss >> contour_id)
      comp.insert(contour_id);

    components.push_back(comp);
  }

  return new Tile(upper, components, parent, is_last);
}

int main() {
  cout << "asdf" << endl;
  cin >> levels;

  // Burn off trailing newline.
  string buffer;
  getline(cin, buffer);

  // Parse in all the tiles.
  Tile **tiles = new Tile*[levels];

  tiles[0] = readTile(NULL, false);

  for (int i = 1; i < levels; i++)
    tiles[i] = readTile(tiles[i-1], (i == levels-1));

  cout << *tiles[levels-1] << endl;
  cout << "is valid: " << tiles[levels-1]->isValid() << endl;
}
