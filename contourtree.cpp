#include <vector>

#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

#include "TilingUtils.h"

int getTets(const char* fn_t, Eigen::MatrixXi &T) {
  FILE* inf_t = fopen(fn_t, "r");
  if (inf_t == NULL) {
    perror("Could not open tet file");
    return -1;
  }

  printf("Getting tets...\n");
  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<std::vector<int> > tets_v;
  while((read = getline(&line, &len, inf_t)) != -1) {
    // Get the tet vertices
    std::vector<int> verts;
    int a, b, c, d;
    int ret;
    if ((ret = sscanf(line, "%d %d %d %d", &a, &b, &c, &d)) != 4) {
      fprintf(stderr, "Error reading line [%s], found %d elem\n", line, ret);
      continue;
    }
    verts.push_back(a);
    verts.push_back(b);
    verts.push_back(c);
    verts.push_back(d);
    tets_v.push_back(verts);
  }
  fclose(inf_t);
  printf("Found %lu tets\n", tets_v.size());

  // Construct the eigen matrix
  T.resize(tets_v.size(), 4);
  for (int i = 0; i < tets_v.size(); ++i) {
    for (int j = 0; j < 4; ++j) {
      T(i,j) = tets_v[i][j];
    }
  }

  if (line) free(line);

  return 0;
}

template<typename T, typename E>
int getOrig(const char* fn, const char* printf_str, E &O, T dummy) {
  FILE* inf = fopen(fn, "r");
  if (!inf) {
    perror("Could not open file");
    return -1;
  }
  printf("Getting lines...\n");

  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<T> orig;
  while((read = getline(&line, &len, inf)) != -1) {
    T val;
    sscanf(line, printf_str, &val);
    orig.push_back(val);
  }
  fclose(inf);

  printf("Found %lu elements\n", orig.size());
  O.resize(orig.size());
  for (int i = 0; i < orig.size(); ++i) {
    O(i) = orig[i];
  }

  return 0;
}
int getHeats(const char* fn_h, Eigen::VectorXd &H) {
  FILE* inf_h = fopen(fn_h, "r");
  if (!inf_h) {
    perror("Could not open heat file");
    return -1;
  }
  printf("Getting heats...\n");

  char * line = NULL;
  size_t len = 0;
  ssize_t read;
  std::vector<double> heats_v;
  while((read = getline(&line, &len, inf_h)) != -1) {
    double val;
    sscanf(line, "%lf", &val);
    heats_v.push_back(val);
  }
  fclose(inf_h);

  printf("Found %lu heat elements\n", heats_v.size());
  H.resize(heats_v.size());
  for (int i = 0; i < heats_v.size(); ++i) {
    H(i) = heats_v[i];
  }

  return 0;
}

int main() {
  Eigen::MatrixXd V;
  Eigen::MatrixXi F, T;
  Eigen::VectorXd H;
  Eigen::VectorXi O;

  const char* fn_off = "debug_tets.off";
  igl::readOFF(fn_off, V,F);
  // Also need to get the tets
  const char* fn_t = "debug_tets.txt";
  int tet = getTets(fn_t, T);
  if (tet < 0) return tet;
  // Also need to get the heat values
  double d;
  int heat = getOrig("debug_heats.txt", "%lf", H, d);
  if (heat < 0) return heat;
  // Also need to get the orig values
  int i;
  int orig = getOrig("debug_orig.txt", "%d", O, i);
  if (orig < 0) return orig;

  igl::viewer::Viewer v;
  v.data.set_mesh(V, F);
  v.launch();

  std::vector<TilingUtils::ConnectedComponent> comps = 
      TilingUtils::allPossibleTiles(V,F,T,O, H);

  // Use dancing links.
}
