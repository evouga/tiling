#include <map>
#include <vector>

#include <igl/remove_duplicates.h>
#include <igl/writeOFF.h>

// Pseudocode copied from http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
//

const double EPS = 1e-6;

struct cmpVectors {
  bool operator()(const std::vector<double> &a, const std::vector<double> &b) {
    if (std::abs(a[0] - b[0]) > EPS) {
      return a[0] - b[0];
    }
    if (std::abs(a[1] - b[1]) > EPS) {
      return a[1] - b[1];
    }
    
    return abs(a[2] - b[2]) > EPS;
  }
};
std::map<std::vector<double>, int, cmpVectors> m_usedPoints;

double m_ConstantSize = 1;

// Note that this will add duplicate points, which need to be removed later.
int addVertex(double x, double y, double z, std::vector<std::vector<double> > &V) {
  // Also make sure the point is on the unit sphere.
  double length = sqrt(x*x + y*y + z*z);

  // Normalize
  std::vector<double> pt = {x/length*m_ConstantSize, y/length*m_ConstantSize, z/length*m_ConstantSize};

  if (m_usedPoints.find(pt) == m_usedPoints.end()) {
    V.push_back(pt);
    m_usedPoints[pt] = V.size() - 1;
    return V.size() - 1;
  } else {
    return m_usedPoints[pt];
  }
}
// Will add a new point to V between points at i and j, returning the index.
int mid(int i, int j, std::vector<std::vector<double> > &V) {
  double x = (V[i][0] + V[j][0]) / 2.0;
  double y = (V[i][1] + V[j][1]) / 2.0;
  double z = (V[i][2] + V[j][2]) / 2.0;

  return addVertex(x, y, z, V);
}

void tile_init(std::vector<std::vector<double> > &V,
               std::vector<std::vector<int> > &F) {

  double t = (1.0 + sqrt(5.0)) / 2.0;

  // Start with corners of three orthogonal rectangles.
  addVertex(-1,  t, 0, V);
  addVertex( 1,  t, 0, V);
  addVertex(-1, -t, 0, V);
  addVertex( 1, -t, 0, V);

  addVertex(0, -1,  t, V);
  addVertex(0,  1,  t, V);
  addVertex(0, -1, -t, V);
  addVertex(0,  1, -t, V);

  addVertex( t, 0, -1, V);
  addVertex( t, 0,  1, V);
  addVertex(-t, 0, -1, V);
  addVertex(-t, 0,  1, V);

  // 5 faces around point 0
  F.push_back({0, 11, 5});
  F.push_back({0, 5, 1});
  F.push_back({0, 1, 7});
  F.push_back({0, 7, 10});
  F.push_back({0, 10, 11});
  // 5 adjacent faces
  F.push_back({1, 5, 9});
  F.push_back({5, 11, 4});
  F.push_back({11, 10, 2});
  F.push_back({10, 7, 6});
  F.push_back({7, 1, 8});

  // 5 faces around point 3
  F.push_back({3, 9, 4});
  F.push_back({3, 4, 2});
  F.push_back({3, 2, 6});
  F.push_back({3, 6, 8});
  F.push_back({3, 8, 9});

  // 5 adjacent faces
  F.push_back({4, 9, 5});
  F.push_back({2, 4, 11});
  F.push_back({6, 2, 10});
  F.push_back({8, 6, 7});
  F.push_back({9, 8, 1});
}

void refine(int k,
            std::vector<std::vector<double> > &V,
            std::vector<std::vector<int> > &F) {
  for (int i = 0; i < k; ++i) {
    printf("  level %d\n", i);
    // New faces.
    std::vector<std::vector<int> > F2;

    for (const std::vector<int> &f : F) {
      int a = mid(f[0], f[1], V);
      int b = mid(f[1], f[2], V);
      int c = mid(f[2], f[0], V);

      
      F2.push_back({f[0], a, c});
      F2.push_back({f[1], b, a});
      F2.push_back({f[2], c, b});
      F2.push_back({a, b, c});
    }

    F = F2;
  }
}

// Will return 1 if v1 < 0 < v2, 2 if v2 < 0 < v1, 0 otherwise.
int sameSign(const std::vector<double> &v1,
             const std::vector<double> &v2) {
  if (v1[0] < 0 && v2[0] > 0) {
    return 1;
  }
  if (v1[0] > 0 && v2[0] < 0) {
    return 2;
  }

  return 0;
}

// Add vertices along the x=0 hemisphere to make the debugging thing more interesting.
void finish(std::vector<std::vector<double> > &V,
            std::vector<std::vector<int> > &F) {

  std::vector<std::vector<int> > new_F;
  for (const auto& f : F) {
    for (int j = 0; j < 3; ++j) {
      int idx0 = f[j];
      int idx1 = f[(j + 1) % 3];
      int idx2 = f[(j + 2) % 3];

      // If they are on different sides of the hemisphere, split them.
      if ( sameSign(V[idx0], V[idx1]) != 0) {
        int newp = mid(idx0, idx1, V);
        // Add two new faces. Clockwise.
        new_F.push_back({idx0, newp, idx2});
        new_F.push_back({newp, idx1, idx2});
      }
    }
  }
  // Add all new faces to F.
  for (const auto& f : new_F) {
    F.push_back(f);
  }
}

void vec_to_Eigen(const std::vector<std::vector<double> > &V,
                  const std::vector<std::vector<int> > &F,
                  Eigen::MatrixXd &eV,
                  Eigen::MatrixXi &eF) {
  // Temp V and F
  Eigen::MatrixXd tV;
  Eigen::MatrixXi tF;
  // Just blindly add all the points in.
  tV.resize(V.size(), 3);
  tF.resize(F.size(), 3);

  for (int i = 0; i < V.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      tV(i, j) = V[i][j];
    }
  }
  for (int i = 0; i < F.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      tF(i, j) = F[i][j];
    }
  }

  //printf("  checking for duplicate points...\n");
  // Map to unique vertices (not used)
  //Eigen::VectorXi toUnique;
  //igl::remove_duplicates(tV, tF, eV, eF, toUnique, 1e-5);
  eV = tV;
  eF = tF;
}


int main(int argc, char* argv[]) {
  if (argc < 4) {
    fprintf(stderr, "Program that will create a sphere from an icosohedral tiling with unit size\n\n");
    fprintf(stderr, "Pseudocode copied from http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html\n\n");
    fprintf(stderr, "usage: %s <level> <size> <output>\n", argv[0]);
    return -1;
  }
  int ref_level = atoi(argv[1]);
  m_ConstantSize = atof(argv[2]);
  const char* out_f = argv[3];

  std::vector<std::vector<double> > V_v;
  std::vector<std::vector<int> > F_v;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  printf("Creating initial icosehedral tiling...\n");
  tile_init(V_v, F_v);
  //vec_to_Eigen(V_v, F_v, V, F);
  //igl::writeOFF("sphere.off", V, F);

  printf("Refining...\n");
  refine(ref_level, V_v, F_v);
  printf("Done! Adding axis points\n");
  // Add points down the x=0 axis.
  finish(V_v, F_v);

  printf("Done! Converting to Eigen and printing out\n");
  vec_to_Eigen(V_v, F_v, V, F);
  igl::writeOFF(out_f, V, F);
  printf("Finshed. See %s for output\n", out_f);
}
