#ifndef TILER_H
#define TILER_H

/* These methods are used to perform the tiling.
 *
 * Currently these are just combinatorics methods.
 * 
 * TODO(bradyz, nclement): implement incremental scoring.
 *
 * */

#include <set>
#include <vector>

namespace Tiler {

/**
 * Each contour is in its own connected component for the first slice.
 *
 * Arguments:
 *  @param number_contours - number of contours on the first slice.
 *
 * Returns:
 *  vector of subsets with a single element.
 */
std::vector<std::set<int> > generateInitialRequirements(int number_contours);

/**
 * Returns the top ids and bottom ids of contours.
 *
 * Arguments:
 *  @param bot_count - number of contours on the bot slice.
 *  @param top_count - number of contours on the top slice.
 *  @param bot_ids - the ids of contours on the bot slice.
 *  @param top_ids - the ids of contours on the top slice.
 *  @param both - the ids of of both contours.
 */
void generateTopAndBottom(int bot_count, int top_count,
                          std::set<int> &bot_ids, std::set<int> &top_ids,
                          std::set<int> &both);

/**
 * Returns all possible subsets.
 *
 * Arguments:
 *  @param total - total number of elements.
 *
 * Returns:
 *  vector of unique subsets.
 */
std::vector<std::set<int> > generatePossibleSubsets(unsigned int total,
                                                    const std::vector<std::set<int> > &previous);

/**
 * Returns if a tiling is valid.
 *
 * This means no loops in the tiling and the previous connected components.
 *
 * Arguments:
 *  @param tile - components used in the tiling.
 *  @param lower - bottom contours to cover.
 *  @param upper - top contours to cover.
 *  @param previous - which components are connected previously.
 *  @param is_last - requires all top components to be connected.
 *
 * Returns:
 *  vector of genus 0 tilings.
 */
bool isValidTiling(const std::vector<std::set<int> > &tile,
                   const std::set<int> &lower,
                   const std::set<int> &upper,
                   const std::vector<std::set<int> > &previous,
                   bool is_last);

/**
 * Populates all valid configurations of connections..
 *
 * Arguments:
 *  @param lower - bottom contours.
 *  @param upper - top contours.
 *  @param all_contours - contours to cover.
 *  @param previous - which components are connected previously.
 *  @param tilings - vector of genus 0 tilings.
 *  @param upper_used - vector of corresponding top pieces.
 *  @param is_last - denotes whether this is the last tile.
 */
void generateTiles(const std::set<int> &lower,
                   const std::set<int> &upper,
                   const std::set<int> &all_contours,
                   const std::vector<std::set<int> > &previous,
                   std::vector<std::vector<std::set<int> > > &tilings,
                   std::vector<std::vector<std::set<int> > > &upper_used,
                   bool is_last);

} // namespace Tiler

#endif
