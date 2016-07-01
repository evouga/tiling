#ifndef TILE_H
#define TILE_H

#include <Eigen/Core>

class Slice;

class Tile
{
public:
	Tile(const Slice &bottom, const Slice &top, bool scale = true);
	~Tile();
	void getOrig(Eigen::MatrixXd &Vtop, Eigen::MatrixXd &Vbot);
	void triangulateSlices(double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
		Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig);

	

private:
	void computeCubeTransformation();
	// Refactored triangualate to add getOrig function
	int addOrig(const Slice &s,
							Eigen::MatrixXd &V, Eigen::MatrixXi &E,
							Eigen::VectorXi &VM, Eigen::VectorXi &EM,
							int offset = 0);

	void triangulateSlice(const Slice &s, 
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double z, double areaBound,
												Eigen::MatrixXd &verts, Eigen::MatrixXi &faces,
												Eigen::VectorXi &orig);

	// Scales a point backward and forward
	void scale(Eigen::Vector2d &pt);
	void unscale(Eigen::Vector2d &pt);

	const static double tilePadding_;

	const Slice &bottom_;
	const Slice &top_;	
	Eigen::Vector2d translate_;
	Eigen::Matrix2d scale_;
  
  bool use_scaling_;
};

#endif
