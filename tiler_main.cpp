#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <cstdlib>

#include "Tiler.h"

using namespace std;
using namespace Tiler;

// Number of contours on top and bottom.
const unsigned int NUMBER_LEVELS = 4;
const unsigned int COUNTS[NUMBER_LEVELS] = {1, 2, 2, 1};

// Holds possible configurations from index 1 and on.
vector<Tile*> generated[NUMBER_LEVELS];

// Generates contour ids.
int total = 0;

int main() {
  vector<Tile*> previous_layer;

  for (int i = 0; i < NUMBER_LEVELS; i++) {
    set<int> upper;
    for (int j = 0; j < COUNTS[i]; j++)
      upper.insert(total++);

    if (i == 0) {
      generated[0].push_back(new Tile(upper));
      continue;
    }

    for (Tile *parent : generated[i-1]) {
      bool last = (i == NUMBER_LEVELS - 1);
      for (Tile *upper_tile : Tiler::generateTiles(upper, parent, last))
        generated[i].push_back(upper_tile);
    }
  }

  for (int i = 0; i < NUMBER_LEVELS; i++)
    cout << "Level " << i << ": " << generated[i].size() << endl;

  cout << endl << "Sample full tile stack top to bottom - " << endl;
  cout << *generated[NUMBER_LEVELS-1][0] << endl;
}
