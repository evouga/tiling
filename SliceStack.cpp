#include <vector>
#include <limits>
#include <algorithm>

#include <igl/colon.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/cotmatrix.h>
#include <igl/jet.h>
#include <igl/min_quad_with_fixed.h>
//#include <igl/point_in_poly.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <igl/triangle/triangulate.h>
#include <igl/viewer/Viewer.h>
#include <igl/writeOFF.h>

#include "glob_defs.h"
#include "offsetSurface.h"
#include "SliceStack.h"
#include "SliceParser.h"
#include "Slice.h"
#include "Tile.h"
#include "viewTetMesh.h"

using namespace std;

SliceStack::SliceStack(const char *baseFilename, const char *objectname) {
  readSlicesFromFolder(baseFilename, objectname, slices_);
  numSlices_ = slices_.size();

  double z = 0.0;
  double xmin = numeric_limits<double>::max();
  double xmax = numeric_limits<double>::min();
  double ymin = numeric_limits<double>::max();
  double ymax = numeric_limits<double>::min();

  for(int i=0; i<numSlices_; i++)
  {
    heights_.push_back(z);
    z += slices_[i]->thickness;
    slices_[i]->xminmax(xmin, xmax);
    slices_[i]->yminmax(ymin, ymax);
  }

  // 3x2
  bbox_.resize(3, 2);
  bbox_.row(0) << xmin, xmax;
  bbox_.row(1) << ymin, ymax;
  bbox_.row(2) << 0, z;
}

SliceStack::~SliceStack() {
  for(int i=0; i<numSlices_; i++)
    delete slices_[i];
}

void SliceStack::triangulateSlice(int start, double areaBound,
                                  Eigen::MatrixXd &botV, Eigen::MatrixXi &botF,
                                  Eigen::MatrixXd &topV, Eigen::MatrixXi &topF,
                                  Eigen::VectorXi &botO, Eigen::VectorXi &topO) {
  assert(start >= 0 && start < slices_.size());

  printf("PRE Number of vertices: ");
  printf("top %lu bot %lu\n", topV.rows(), botV.rows());

  Tile t(*slices_[start], *slices_[start+1], bbox_);
  t.triangulateSlices(areaBound, botV, botF, topV, topF, botO, topO);

  printf("Number of vertices: top %lu bot %lu\n", topV.rows(), botV.rows());

  // Then flip normals of bottom slice
  flipNormal(botF);
}

int SliceStack::getSizeAt(int i) {
  if (i > slices_.size())
    return -1;
  return slices_[i]->contours.size();
}

bool customSortByX(const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[0] < b[0];
}

bool customSortByY(const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[1] < b[1];
}

bool customSortByZ(const Eigen::Vector3d a, const Eigen::Vector3d b) {
  return a[2] < b[2];
}

//                      left             right              back          front
// constantCoord => 0: (-0.5, _, _), 1: (0.5, _, _), 2: (_, -0.5, _), 3: (_, 0.5, _)
// also include the minimum and maximum values in the non-constant coord (mn, mx)
// and the minimum and maximum values for this coord (o_mn, o_mx)
void SliceStack::triangulateSide(int constantCoord,
                                 double mn, double mx,
                                 double o_mn, double o_mx,
                                 vector<Eigen::Vector3d> &verts,
                                 Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

  sort(verts.begin(), verts.end(), customSortByZ);

  vector<Eigen::Vector3d> topV;
  vector<Eigen::Vector3d> botV;

  // Compute the mid z-coord
  double mid = 0;
  for (const auto & v : verts) {
    mid += v(2);
  }
  mid /= verts.size();

  // Fill the botV and topV vectors
  for (int i = 0; i < verts.size(); ++i) {
    if (verts[i][2] < mid)
      botV.push_back(verts[i]);
    else
      topV.push_back(verts[i]);
  }

  // Set the z-scaling as the distance between the top and the bottom vert
  double z_spacing = verts.back()(2) - verts.front()(2);

  // If the constant coord is y, sort by x direction.
  if (constantCoord >= 2) {
    sort(botV.begin(), botV.end(), customSortByX);
    sort(topV.rbegin(), topV.rend(), customSortByX); // backwards
  }
  // otherwise, sort by y direction.
  else {
    sort(botV.begin(), botV.end(), customSortByY);
    sort(topV.rbegin(), topV.rend(), customSortByY); // backwards
  }

  Eigen::MatrixXd inputV(verts.size() + GLOBAL::EXTRA*2, 2);
  Eigen::MatrixXi inputE(verts.size() + GLOBAL::EXTRA*2, 2);
  Eigen::MatrixXd H(0, 2);

  int offset = 0;

  // Add bottom points, from left to right
  for (int i = 0; i < botV.size(); ++i) {
    if (constantCoord < 2) {
      inputV.row(offset++) << botV[i][1], botV[i][2];
    } else {
      inputV.row(offset++) << botV[i][0], botV[i][2];
    }
  }

  // Add right points
  // Get the vector that extends from the top corner to the bottom corner.
  Eigen::Vector3d bot2Top = topV.front() - botV.back(); // top is sorted backwards
  for (int i = 1; i <= GLOBAL::EXTRA; ++i) {
    auto shift_by = botV.back() + bot2Top * i / (GLOBAL::EXTRA + 1);
    if (constantCoord >= 2) {  //constant is y
      inputV.row(offset++) << shift_by(0), shift_by(2);
    } else {  // constant is x
      inputV.row(offset++) << shift_by(1), shift_by(2);
    }
    //inputV.row(offset++) = 
        //Eigen::Vector2d(mx, botV[0][2] + z_spacing / GLOBAL::EXTRA * i);
  }

  // Add top points, from right to left
  for (int i = 0; i < topV.size(); ++i) {
    if (constantCoord >= 2) {
      inputV.row(offset++) << topV[i][0], topV[i][2];
    } else {
      inputV.row(offset++) << topV[i][1], topV[i][2];
    }
  }

  // Add left points
  Eigen::Vector3d top2Bot = botV.front() - topV.back();
  for (int i = 1; i <= GLOBAL::EXTRA; ++i) {
    auto shift_by = topV.back() + top2Bot * i / (GLOBAL::EXTRA + 1);
    if (constantCoord >= 2) {  //constant is y
      inputV.row(offset++) << shift_by(0), shift_by(2);
    } else {  // constant is x
      inputV.row(offset++) << shift_by(1), shift_by(2);
    }
    //inputV.row(offset++) =
        //Eigen::Vector2d(mn, topV[0][2] - z_spacing / GLOBAL::EXTRA * i);
  }

  for (int i = 0; i < inputV.rows(); ++i) {
    inputE.row(i) << i, (i + 1) % (inputV.rows());
  }

  Eigen::MatrixXd tmpV;

	cout << "Delaunay triangulating mesh with " << inputV.rows() << " verts" << endl;
  stringstream params;
  params << "DYa" << triangle_max_area;
  params << "q";
  cout << "parameters are " << params.str() << endl;
	igl::triangle::triangulate(inputV,inputE,H,params.str().c_str(),tmpV,F);	
	cout << "done" << endl;

  V.resize(tmpV.rows(), 3);

  for (int i = 0; i < tmpV.rows(); ++i) {
    if (constantCoord == 0)
      V.row(i) << o_mn, tmpV(i, 0), tmpV(i, 1);
    else if (constantCoord == 1)
      V.row(i) << o_mx, tmpV(i, 0), tmpV(i, 1);
    else if (constantCoord == 2)
      V.row(i) << tmpV(i, 0), o_mn, tmpV(i, 1);
    else if (constantCoord == 3)
      V.row(i) << tmpV(i, 0), o_mx, tmpV(i, 1);
  }
}

void SliceStack::flipNormal(Eigen::MatrixXi &f) {
	for (int i = 0; i < f.rows(); ++i) {
		int temp = f(i,0);
		f(i,0) = f(i,2);
		f(i,2) = temp;
	}
}

void SliceStack::tetrahedralizeSlice (
        const Eigen::MatrixXd &botV, const Eigen::MatrixXi &botF,
        const Eigen::MatrixXd &topV, const Eigen::MatrixXi &topF,
        const Eigen::VectorXi &botorig, const Eigen::VectorXi &toporig,
        Eigen::MatrixXd &TV, Eigen::MatrixXi &TT,
        Eigen::MatrixXi &TF, Eigen::VectorXi &TO) {
  auto botMin = botV.colwise().minCoeff();
  auto botMax = botV.colwise().maxCoeff();
  auto topMin = topV.colwise().minCoeff();
  auto topMax = topV.colwise().maxCoeff();

  double xmin = std::min(botMin(0), topMin(0)),
         ymin = std::min(botMin(1), topMin(1)),
         zmin = std::min(botMin(2), topMin(2)),
         xmax = std::max(botMax(0), topMax(0)),
         ymax = std::max(botMax(1), topMax(1)),
         zmax = std::max(botMax(2), topMax(2));

  printf("bounding box is %f,%f %f,%f %f,%f\n", xmin,xmax, ymin,ymax, zmin,zmax);

  vector<Eigen::Vector3d> leftV;
  vector<Eigen::Vector3d> rightV;
  vector<Eigen::Vector3d> frontV;
  vector<Eigen::Vector3d> backV;

  for (int i = 0; i < botV.rows(); ++i) {
    double x = botV(i, 0),
           y = botV(i, 1);
    // x coordinate touches
    if (x <= botMin(0) + GLOBAL::EPS)
      leftV.push_back(botV.row(i));
    if (x >= botMax(0) - GLOBAL::EPS)
      rightV.push_back(botV.row(i));
    // y coordinate touches
    if (y <= botMin(1) + GLOBAL::EPS)
      frontV.push_back(botV.row(i));
    if (y >= botMax(1) - GLOBAL::EPS)
      backV.push_back(botV.row(i));
  }

  for (int i = 0; i < topV.rows(); ++i) {
    double x = topV(i, 0),
           y = topV(i, 1);
    // x coordinate touches
    if (x <= topMin(0) + GLOBAL::EPS)
      leftV.push_back(topV.row(i));
    if (x >= topMax(0) - GLOBAL::EPS)
      rightV.push_back(topV.row(i));
    // y coordinate touches
    if (y <= topMin(1) + GLOBAL::EPS)
      frontV.push_back(topV.row(i));
    if (y >= topMax(1) - GLOBAL::EPS)
      backV.push_back(topV.row(i));
  }

  Eigen::MatrixXd vall;
  vall.resize(leftV.size() + rightV.size() + frontV.size() + backV.size(), 3);
  Eigen::VectorXi typeall;
  typeall.resize(vall.rows());
  for (int i = 0; i < vall.rows(); ++i) {
    if (i < leftV.size()) {
      typeall(i) = 0;
      vall.row(i) = leftV[i];
    } else if (i < leftV.size() + rightV.size()) {
      typeall(i) = 1;
      vall.row(i) = rightV[i - leftV.size()];
    } else if (i < leftV.size() + rightV.size() + frontV.size()) {
      typeall(i) = 2;
      vall.row(i) = frontV[i - rightV.size() - leftV.size()];
    } else {
      typeall(i) = 3;
      vall.row(i) = backV[i - leftV.size() - rightV.size() - frontV.size()];
    }
  }
  Eigen::MatrixXd colsall;
	igl::jet(typeall, true, colsall);
  igl::viewer::Viewer v;
  v.data.set_points(vall, colsall);
  v.launch();

  Eigen::MatrixXd leftTriV;
  Eigen::MatrixXi leftTriF;

  Eigen::MatrixXd rightTriV;
  Eigen::MatrixXi rightTriF;

  Eigen::MatrixXd backTriV;
  Eigen::MatrixXi backTriF;

  Eigen::MatrixXd frontTriV;
  Eigen::MatrixXi frontTriF;

  triangulateSide(0, ymin,ymax, xmin,xmax, leftV, leftTriV, leftTriF);
	flipNormal(leftTriF);
  v.data.set_mesh(leftTriV, leftTriF);
  v.launch();
  triangulateSide(2, xmin,xmax, ymin,ymax, frontV, frontTriV, frontTriF);
  v.data.clear();
  v.data.set_mesh(frontTriV, frontTriF);
  v.launch();
  triangulateSide(1, ymin,ymax, xmin,xmax, rightV, rightTriV, rightTriF);
  v.data.clear();
  v.data.set_mesh(rightTriV, rightTriF);
  v.launch();
  triangulateSide(3, xmin,xmax, ymin,ymax, backV, backTriV, backTriF);
	flipNormal(backTriF);
  v.data.clear();
  v.data.set_mesh(backTriV, backTriF);
  v.launch();

	// Can't count points duplicate times

  int totalVertices = topV.rows() + botV.rows() + 
    leftTriV.rows() + rightTriV.rows() + backTriV.rows() + frontTriV.rows();

  int totalFaces = topF.rows() + botF.rows() + 
    leftTriF.rows() + rightTriF.rows() + backTriF.rows() + frontTriF.rows();

  Eigen::MatrixXd V_rep(totalVertices, 3);
	Eigen::VectorXi M_rep(totalVertices);

  int offset = 0;

  // Add the vertices
  for (int i = 0; i < botV.rows(); ++i) {
		M_rep(offset) = botorig(i);
    V_rep.row(offset++) = botV.row(i);
  }

  for (int i = 0; i < topV.rows(); ++i) {
		M_rep(offset) = toporig(i);
    V_rep.row(offset++) = topV.row(i);
  }

  for (int i = 0; i < leftTriV.rows(); ++i) {
		// Not original
		M_rep(offset) = GLOBAL::nonoriginal_marker;
    V_rep.row(offset++) = leftTriV.row(i);
  }

  for (int i = 0; i < rightTriV.rows(); ++i) {
		// Not original
		M_rep(offset) = GLOBAL::nonoriginal_marker;
    V_rep.row(offset++) = rightTriV.row(i);
  }

  for (int i = 0; i < frontTriV.rows(); ++i) {
		// Not original
		M_rep(offset) = GLOBAL::nonoriginal_marker;
    V_rep.row(offset++) = frontTriV.row(i);
  }

  for (int i = 0; i < backTriV.rows(); ++i) {
		// Not original
		M_rep(offset) = GLOBAL::nonoriginal_marker;
    V_rep.row(offset++) = backTriV.row(i);
  }

	// Get all the unique vertices
	Eigen::MatrixXd V;
	// contains mapping from unique to all indices
	// Size: #V
	Eigen::VectorXi unique_to_all;
	// contains mapping from all indices to unique
	// Size: #V_rep
	Eigen::VectorXi all_to_unique;
	igl::unique_rows(V_rep, V,unique_to_all,all_to_unique);
	//printf("Size of unique is now: %ld vs %ld\n", V_rep.rows(), V.rows());
	//printf("Sizes are %ld,%ld\n", unique_to_all.rows(), all_to_unique.rows());

	// Get unique markers for M
	Eigen::VectorXi M(V.rows());
	for (int i = 0; i < M.rows(); ++i) {
		//printf("Changing index %d to %d\n", i, unique_to_all(i));
		M(i) = M_rep(unique_to_all(i));
	}
  
	// Add the faces
  Eigen::MatrixXi F(totalFaces, 3);
  Eigen::Vector3i vertexOffset(0, 0, 0);
  offset = 0;

  for (int i = 0; i < botF.rows(); ++i) {
    F.row(offset++) = botF.row(i);
  }

  vertexOffset += Eigen::Vector3i(botV.rows(), botV.rows(), botV.rows());

  for (int i = 0; i < topF.rows(); ++i) {
    F.row(offset++) = vertexOffset + Eigen::Vector3i(topF.row(i));
  }

  vertexOffset += Eigen::Vector3i(topV.rows(), topV.rows(), topV.rows());

  for (int i = 0; i < leftTriF.rows(); ++i) {
    F.row(offset++) = vertexOffset + Eigen::Vector3i(leftTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(leftTriV.rows(), leftTriV.rows(), leftTriV.rows());

  for (int i = 0; i < rightTriF.rows(); ++i) {
    F.row(offset++) = vertexOffset + Eigen::Vector3i(rightTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(rightTriV.rows(), rightTriV.rows(), rightTriV.rows());

  for (int i = 0; i < frontTriF.rows(); ++i) {
    F.row(offset++) = vertexOffset + Eigen::Vector3i(frontTriF.row(i));
  }

  vertexOffset += Eigen::Vector3i(frontTriV.rows(), frontTriV.rows(), frontTriV.rows());

  for (int i = 0; i < backTriF.rows(); ++i) {
    F.row(offset++) = vertexOffset + Eigen::Vector3i(backTriF.row(i));
  }

	// Make sure the faces point to the correct (unique) points
	for (int i = 0; i < F.rows(); ++i) {
		for (int j = 0; j < F.row(i).cols(); ++j) {
			F(i,j) = all_to_unique(F(i,j));
		}
	}
	

  //igl::writeOFF("foo.off", V, F);

	Eigen::VectorXi FM(F.rows());
	for (int i = 0; i < F.rows(); ++i) {
		if (M(F(i,0)) == GLOBAL::original_marker && 
        M(F(i,1)) == GLOBAL::original_marker && 
        M(F(i,2)) == GLOBAL::original_marker) {
			FM(i) = GLOBAL::original_marker;
		} else {
			FM(i) = GLOBAL::nonoriginal_marker;
		}
	}
  // Tetrahedralized interior

  //igl::viewer::Viewer v;
  v.data.clear();
  v.data.set_mesh(V, F);
  v.launch();

	// TV will have the tetrahedralized vertices;
	// TT will have the "" tet indices (#V x 4)
	// TF will have the "" face indices (#V x 3)
	// TO will have the "" vertex markers
  stringstream params;
  params << "pq" << tetgen_max_rad_ratio;
  params << "Y";
  igl::copyleft::tetgen::tetrahedralize(V,F,M,FM, params.str().c_str(), TV,TT,TF,TO);
  igl::writeOFF("foo_tet.off", TV, TF);
  printf("Tetrahedralize done\n");
  printf("Number of faces before:%lu and after:%lu\n", F.rows(), TF.rows());
}

void SliceStack::computeLaplace(int slice_no,
																const Eigen::MatrixXd &TV,
																const Eigen::MatrixXi &TT,
																const Eigen::MatrixXi &TF,
																const Eigen::VectorXi &TO,
                                Eigen::VectorXd &Z) {
  bool laplace_DEBUG = true;

	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	Eigen::IOFormat LongFmt(10, 0, ", ", "\n", "[", "]");
	Eigen::IOFormat RFmt(4, 0, ", ", ", ", "", "", "(", ")");
	assert(TO.rows() == TV.rows());

  if (laplace_DEBUG) printf("Setting constant values...");
	Eigen::VectorXi known;
	Eigen::VectorXd known_c;
	std::vector<int> known_v;
	std::vector<double> known_c_v;
  auto mx = TV.colwise().maxCoeff();
  auto mn = TV.colwise().minCoeff();

	for (int i = 0; i < TO.rows(); ++i) {
		/*
		if (TV(i,2) == -GLOBAL::z_lim) {
			known_v.push_back(i);
			known_c_v.push_back(1);
		} else if (TV(i,2) == GLOBAL::z_lim) {
			known_v.push_back(i);
			known_c_v.push_back(0);
		}
		*/
		if (TO(i) == GLOBAL::original_marker) {
			known_v.push_back(i);
			known_c_v.push_back(GLOBAL::inside_temp);
		} else if (TV(i,2) == mx(2) || TV(i,2) == mn(2)) {
			known_v.push_back(i);
			known_c_v.push_back(GLOBAL::outside_temp);
		}
	}
	if (laplace_DEBUG) 
    printf("done! Number of known values is %lu/%lu\n",
           known_v.size(), TV.rows());
	known.resize(known_v.size());
	known_c.resize(known_v.size());
	for (int i = 0; i < known_c.size(); ++i) {
		known(i) = known_v[i];
		known_c(i) = known_c_v[i];
	}
	
  if (laplace_DEBUG) 
    printf("Constructing Laplacian...");
	// Construct Laplacian
	// Dense matrix first
	//
	Eigen::SparseMatrix<double> L(TV.rows(), TV.rows());
	// Set non-diag elements to 1 if connected, 0 otherwise
	// Use the tets instead of the faces
	for (int i = 0; i < TT.rows(); ++i) {
		L.coeffRef(TT(i,0), TT(i,1)) = -1; L.coeffRef(TT(i,1), TT(i,0)) = -1;
		L.coeffRef(TT(i,1), TT(i,2)) = -1; L.coeffRef(TT(i,2), TT(i,1)) = -1;
		L.coeffRef(TT(i,2), TT(i,3)) = -1; L.coeffRef(TT(i,3), TT(i,2)) = -1;
		L.coeffRef(TT(i,3), TT(i,0)) = -1; L.coeffRef(TT(i,0), TT(i,3)) = -1;
	}

	// Set diag elements to valence of entry
	for (int i = 0; i < TV.rows(); ++i) {
		L.coeffRef(i,i) = -L.row(i).sum();
	}
	//Eigen::SparseMatrix<double> L = L_d.sparseView();//, L_known, L_others;
  if (laplace_DEBUG)
    printf("done! Number non-zeros is %lu\n", L.nonZeros());
	
  if (laplace_DEBUG)
    printf("Solving energy constraints...");

	// Solve energy constraints.
	igl::min_quad_with_fixed_data<double> mqwf;
	// Linear term is 0
	Eigen::VectorXd B = Eigen::VectorXd::Zero(TV.rows(), 1);
	// Empty Constraints
	Eigen::VectorXd Beq;
	Eigen::SparseMatrix<double> Aeq;
	bool success = igl::min_quad_with_fixed_precompute(L, known, Aeq, false, mqwf);
	if(!success)
		fprintf(stderr,"ERROR: fixed_precompute didn't work!\n");
	igl::min_quad_with_fixed_solve(mqwf,B,known_c,Beq, Z);

  if (laplace_DEBUG)
    printf("done!\n");

  // Here's how to view it all:
  /*
	// Pseudo-color based on solution
	Eigen::MatrixXd C;
	igl::jet(Z, true, C);
  printf("Min Z is %lf and max Z is %lf\n", Z.minCoeff(), Z.maxCoeff());

  TetMeshViewer::viewTetMesh(TV, TT, TF, Z, true);

  // Plot the mesh with pseudocolors
  // igl::viewer::Viewer viewer;
  // viewer.data.set_mesh(TV, TF);
	// viewer.data.set_face_based(true);
  // viewer.core.show_lines = false;
  // viewer.data.set_colors(C);
  // viewer.launch();

  igl::writeOFF("triangulation.off", TV, TF);
  TetMeshViewer::viewOffsetSurface(TV, TF, TT, Z);
  */
}
