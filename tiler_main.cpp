#include <map>
#include <set>
#include <vector>
#include <iostream>

#include "Tiler.h"

using namespace std;

// Number of contours on top and bottom.
const unsigned int COUNTS[] = {1, 3, 2};
const unsigned int NUMBER_LEVELS = 3;

// Holds possible configurations from index 1 and on.
vector<vector<set<int> > > result[NUMBER_LEVELS];

// Holds what components are connected for the corresponding result.
vector<vector<set<int> > > previous_connected[NUMBER_LEVELS];

map<int, int> generated_from[NUMBER_LEVELS];

int main() {
  // Seed the initial previous connected components.
  previous_connected[0].push_back(Tiler::generateInitialRequirements(COUNTS[0]));

  for (int i = 1; i < NUMBER_LEVELS; i++) {
    int lower_count = COUNTS[i-1];
    int upper_count = COUNTS[i];
    int total = lower_count + upper_count;

    set<int> lower;
    set<int> upper;
    set<int> all_contours;
    Tiler::generateTopAndBottom(lower_count, upper_count,
                                lower, upper, all_contours);

    int tiles_generated = 0;
    for (int j = 0; j < previous_connected[i-1].size(); j++) {
      const vector<set<int> > &connected = previous_connected[i-1][j];
      // A tiling is a set of connected components that cover all contours.
      vector<vector<set<int> > > tilings;
      vector<vector<set<int> > > tilings_upper;
      Tiler::generateTiles(lower, upper, all_contours, connected,
                           tilings, tilings_upper, (i == (NUMBER_LEVELS - 1)));

      for (int k = 0; k < tilings.size(); k++) {
        const vector<set<int> > &tile = tilings[k];
        const vector<set<int> > &tile_upper = tilings_upper[k];

        generated_from[i][tiles_generated] = j;
        result[i].push_back(tile);
        previous_connected[i].push_back(tile_upper);
        tiles_generated++;
      }
    }
  }

  int level = NUMBER_LEVELS - 1;
  for (int i = result[level].size() - 10; i < result[level].size(); i++) {
    int current = i;
    int current_level = level;
    while (current_level > 0) {
      cout << "level: " << current_level << endl;
      const vector<set<int> > &tile = result[current_level][current];
      for (const set<int> &component : tile) {
        for (int x : component)
          cout << x << " ";
        cout << endl;
      }
      current = generated_from[current_level][current];
      current_level--;
    }
    cout << endl;
  }

  cout << result[level].size() << " total number of correct meshes. " << endl;
}
