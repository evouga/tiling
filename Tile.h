#ifndef TILE_H
#define TILE_H

#include <Eigen/Core>

class Slice;

class Tile {
public:
	Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox);
	~Tile();

	void getOrig(Eigen::MatrixXd &Vtop, Eigen::MatrixXd &Vbot);
	void triangulateSlices(double areaBound,
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces,
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
		Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig);

private:
	// Refactored triangualate to add getOrig function
	int addOrig(const Slice &s,
							Eigen::MatrixXd &V, Eigen::MatrixXi &E,
							Eigen::VectorXi &VM, Eigen::VectorXi &EM,
							Eigen::MatrixXd &lims, int offset = 0);

	void triangulateSlice(const Slice &s,
                        double z, double areaBound,
												Eigen::MatrixXd &verts, Eigen::MatrixXi &faces,
												Eigen::VectorXi &orig);

	const static double tilePadding_;

	const Slice &bottom_;
	const Slice &top_;
  Eigen::MatrixXd bbox_;
};

#endif
