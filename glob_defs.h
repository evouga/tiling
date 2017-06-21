#ifndef GLOB_DEFS_H

// Avoid the use of magic numbers and multiple defines.
// Put all the globally-used numbers in here.
namespace GLOBAL {
  // Print/viewer debug mode.
  const bool DEBUG = true;

  const float padding_perc = 0.05;

  // Values defined on the min- and max-z-axis extremes. Inside temp is the
  // temperature inside the contour, and outside temp is the value outside.
  // outside_temp must be greater than inside_temp if the offset surfacing
  // is going to work.
  const float inside_temp = 0.f;
  const float outside_temp = 1.f;

  // This is the upper and lower z-axis amounts
  const float z_lim = 0.5;

  // Used in offsetSurface to define the 'temperature' outside the volume.
  // Must be anything greater than outside_temp.
  const float highest_temp = 100.f;

  // Original marker, used by tetgen and triangulate.
  const int original_marker = 2;
  // nonoriginal_marker is used by tetgen for defining non-contour faces
  const int nonoriginal_marker = 1;
  const int extended_vertices_contour = -1;
  const int nonmanifold_marker = 1000;

  // Used to compare doubles.
  const double EPS = 1e-6;

  // Number of extra points to add to sides of the unit cube.
  // TODO: relate this to thickness.
  const int EXTRA = 20;

  // Labels for triangulating a rhombus.
  const int LEFT = 0;
  const int RIGHT = 1;
  const int FRONT = 2;
  const int BACK = 3;
  
  // Variable used by triangle to determine the maximum area of each triangle.
  // See -a flag. Lower value creates smaller triangles. With "hack" in place
  // (see Tile.cpp), a value of 0.04 isn't terrible.
  //const float triangle_max_area = 0.04;
  const float TRI_AREA = 0.004;
  // Variable used by tetgen to determine the maximum allowed *ratio* between
  // the radius and the area. See -q flag. With "hack" in place (see Tile.cpp),
  // a value of 1.4 isn't terrible.
  //const float tetgen_max_rad_ratio = 1.4;
  const float TET_RATIO = 1.1;
}

#endif
