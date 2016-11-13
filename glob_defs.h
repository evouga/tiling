#ifndef GLOB_DEFS_H

// Avoid the use of magic numbers and multiple defines. Put all the globally-
// used numbers in here.
namespace GLOBAL {
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

  // Used to compare doubles.
  const double EPS = 1e-5;

  // Number of extra points to add to sides of the unit cube.
  // TODO: relate this to thickness.
  const int EXTRA = 20;

  // Labels for triangulating a rhombus
  const int LEFT = 0;
  const int RIGHT = 1;
  const int FRONT = 2;
  const int BACK = 3;
}

#endif
