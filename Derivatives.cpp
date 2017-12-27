#include "Derivatives.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <Eigen/Geometry>

#include <igl/cotmatrix_entries.h>

const double PI = 3.141592653589793238462643383279502884L;
typedef Eigen::Matrix<double, -1, -1, Eigen::RowMajor> RowMatrixXd;
typedef Eigen::Matrix<int, -1, -1, Eigen::RowMajor> RowMatrixXi;

// Functions that are not needed to be exported globally.
namespace {
void accumulateStar(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::MatrixXi &E,
                           Eigen::VectorXd *oneStar /* optional */,
                           std::vector<Eigen::Triplet<double> > *dOneStar /* optional */,
                           int edgeid, int faceid, int facevert) {
  Eigen::Vector3d qi = V.row(F(faceid, facevert));
  Eigen::Vector3d qj = V.row(F(faceid, (facevert+1)%3));
  Eigen::Vector3d qk = V.row(F(faceid, (facevert+2)%3));
  Eigen::Vector3d n = (qj - qk).cross(qi - qk);
  double nnorm = n.norm();
  double cijk = (qj - qk).dot(qi - qk) / nnorm;
  if (oneStar)
    (*oneStar)[edgeid] += 0.5*cijk;

  if (dOneStar)
  {
    Eigen::Vector3d di = 1.0 / nnorm * (qj - qk) - (qj - qk).dot(qi - qk) / nnorm / nnorm / nnorm * n.cross(qj - qk);
    Eigen::Vector3d dj = 1.0 / nnorm * (qi - qk) - (qj - qk).dot(qi - qk) / nnorm / nnorm / nnorm * (qi - qk).cross(n);
    Eigen::Vector3d dk = -di - dj;
    for (int i = 0; i < 3; i++)
    {
      dOneStar->push_back(Eigen::Triplet<double>(edgeid, 3 * F(faceid, facevert) + i, 0.5*di[i]));
      dOneStar->push_back(Eigen::Triplet<double>(edgeid, 3 * F(faceid, (facevert + 1) % 3) + i, 0.5*dj[i]));
      dOneStar->push_back(Eigen::Triplet<double>(edgeid, 3 * F(faceid, (facevert + 2) % 3) + i, 0.5*dk[i]));
    }
  }
}

} // namespace

namespace derivatives {
/**
 * Decompose L into D^T\star D.
 *
 * @param V - vertices of mesh
 * @param F - faces of mesh
 */
void decompose_L(const RowMatrixXd &V, const RowMatrixXi &F,
                 Eigen::SparseMatrix<double> &D, Eigen::SparseMatrix<double> &star) {
  // Construct D
  D.resize(F.rows() * 3, V.rows());
  //                                 # non-zeros per vertex
  D.reserve(Eigen::VectorXi::Constant(2, V.rows()));
  //D.reserve(F.rows()*3 * 2);
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int from = F(i, j);
      int to = F(i, (j + 1) % 3);
      D.insert(i*3 + j, from) = -1;
      D.insert(i*3 + j, to) = 1;
    }
  }

  // Construct star from cot_entries.
  Eigen::MatrixXd cot_entries;
  // cot_entries is #Fx3 where each row corresponds to edge:
  //   ( [1,2], [2,0], [0,1] )
  igl::cotmatrix_entries(V,F, cot_entries);

  star.resize(F.rows() * 3, F.rows() * 3);
  star.reserve(Eigen::VectorXi::Constant(F.rows() * 3, 1));
  for (int i = 0; i < F.rows(); ++i) {
    for (int j = 0; j < 3; ++j) {
      int e = i * 3 + j;
      // Need the negative value for this to work.
      star.coeffRef(e, e) -= cot_entries(i, (j + 2) % 3);
    }
  }

  //fprintf(stderr, "D is %ld,%ld and star is %ld,%ld\nFinished with matrix\n",
  //        D.rows(), D.cols(), star.rows(), star.cols());
}

// Will build edges needed for these methods.
void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E) {
  std::map<std::pair<int, int>,
      Eigen::Vector4i, std::less<std::pair<int, int> >, 
      Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >> edgemap;

  int nfaces = F.rows();
  for (int i = 0; i < nfaces; i++) {
    for (int j = 0; j < 3; j++) {
      int idx1 = F(i, j);
      int idx2 = F(i, (j + 1) % 3);
      int slot = 0;

      if (idx1 > idx2) {
        std::swap(idx1, idx2);
        slot = 1;
      }
      std::map<std::pair<int, int>,
          Eigen::Vector4i, std::less<std::pair<int, int> >,
          Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>::iterator it =
              edgemap.find(std::pair<int, int>(idx1, idx2));
      if (it == edgemap.end()) {
        Eigen::Vector4i newedge;
        newedge[0] = idx1;
        newedge[1] = idx2;
        newedge[2] = newedge[3] = -1;
        newedge[2 + slot] = i;
        edgemap[std::make_pair(idx1, idx2)] = newedge;
      }
      else
      {
        edgemap[std::make_pair(idx1, idx2)][2 + slot] = i;
      }
    }
  }

  int nedges = edgemap.size();
  E.resize(nedges, 4);
  int idx = 0;
  for (std::map<std::pair<int, int>, Eigen::Vector4i>::iterator it = edgemap.begin();
       it != edgemap.end();
       ++it) {
    E.row(idx) = it->second.transpose();
    idx++;
  }
}

void hodgeOneStar(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                  const Eigen::MatrixXi &E,
                  Eigen::VectorXd *oneStar, /* optional */
                  std::vector<Eigen::Triplet<double> > *dOneStar /* optional */) {
  int nfaces = F.rows();
  int nedges = E.rows();
  if (oneStar) {
    oneStar->resize(nedges);
    oneStar->setZero();
  }
  if (dOneStar)
    dOneStar->clear();

  for (int i = 0; i < nedges; i++) {
    int leftface = E(i, 2);
    if (leftface != -1) {
      int idx = -1;
      for (int j = 0; j < 3; j++) {
        if (F(leftface, j) == E(i, 0))
          idx = j;
      }
      assert(idx != -1);
      accumulateStar(V, F, E, oneStar, dOneStar, i, leftface, idx);
    }
    int rightface = E(i, 3);
    if (rightface != -1) {
      int idx = -1;
      for (int j = 0; j < 3; j++) {
        if (F(rightface, j) == E(i, 1))
          idx = j;
      }
      assert(idx != -1);
      accumulateStar(V, F, E, oneStar, dOneStar, i, rightface, idx);
    }
  }
}

void differentialMatrix(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E,
                        std::vector<Eigen::Triplet<double> > &dCoeffs) {
  dCoeffs.clear();
  int nedges = E.rows();
  for (int i = 0; i < nedges; i++) {
    dCoeffs.push_back(Eigen::Triplet<double>(i, E(i, 0), 1.0));
    dCoeffs.push_back(Eigen::Triplet<double>(i, E(i, 1), -1.0));
  }
}

double geodesicCurvature(const RowMatrixXd &V, const RowMatrixXi &F,
                         const RowMatrixXi &E,
                         const std::vector<int> *border_vertices,
                         RowMatrixXd *dGeodesicCurvature /* optional */) {
  if (dGeodesicCurvature) {
    dGeodesicCurvature->resize(V.rows(), 3);
    dGeodesicCurvature->setZero();
  }

  int nedges = E.rows();
  std::set<int> bdryverts;
  double tottheta = 0;

  if (border_vertices) {
    bdryverts.insert(border_vertices->begin(), border_vertices->end());
  } else {
    for (int i = 0; i < nedges; i++) {
      int faceid = -1;
      if (E(i, 2) == -1 || E(i,3) == -1) {
        bdryverts.insert(E(i, 0));
        bdryverts.insert(E(i, 1));            
      }
    }
  }

  int nfaces = F.rows();

  for (int i = 0; i < nfaces; i++) {
    for (int j = 0; j < 3; j++) {
      if (bdryverts.count(F(i, j))) {
        int nextj = (j + 1) % 3;
        int prevj = (j + 2) % 3;
        Eigen::Vector3d p = V.row(F(i, j));
        Eigen::Vector3d nextp = V.row(F(i, nextj));
        Eigen::Vector3d prevp = V.row(F(i, prevj));
        Eigen::Vector3d facenormal = (nextp - p).cross(prevp - p);
        double norm = facenormal.norm();
        double theta = 2.0 * atan2(norm, (nextp - p).norm() * (prevp - p).norm() + (nextp - p).dot(prevp - p));
        tottheta += theta;
        if (dGeodesicCurvature) {
          Eigen::Vector3d z = facenormal / norm;
          Eigen::Vector3d dnext = (nextp - p).cross(z) / (nextp - p).dot(nextp - p);
          Eigen::Vector3d dprev = -(prevp - p).cross(z) / (prevp - p).dot(prevp - p);
          dGeodesicCurvature->row(F(i, nextj)) += -dnext;
          dGeodesicCurvature->row(F(i, prevj)) += -dprev;
          dGeodesicCurvature->row(F(i, j)) += (dnext + dprev);
        }
      }
    }
  }
  return bdryverts.size()*PI - tottheta;
}
// Helper.
double geodesicCurvature(const RowMatrixXd &V, const RowMatrixXi &F,
                         const std::vector<int> *border_vertices,
                         RowMatrixXd *dGeodesicCurvature /* optional */) {
  Eigen::MatrixXi E;
  buildEdges(F, E);
  return geodesicCurvature(V, F, E, border_vertices, dGeodesicCurvature);
}


} // namespace derivatives


