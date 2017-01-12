#ifndef TILE_H
#define TILE_H

#include <Eigen/Core>

class Slice;

class Tile {

public:
	Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox);
	Tile(const Slice &bottom, const Slice &top, const Eigen::MatrixXd &bbox,
       const std::vector<int> &allowed_bot, const std::vector<int> &allowed_top);
	~Tile();

	void triangulateSlices(double areaBound,
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces,
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces,
		Eigen::VectorXi &bot_orig, Eigen::VectorXi &top_orig);

private:
	void addOrig(const Slice &s,
              std::vector<Eigen::RowVector2d> &V,
              std::vector<Eigen::RowVector2i> &E,
              std::vector<int> &VM, std::vector<int> &EM,
							Eigen::MatrixXd &lims,
              const std::vector<int> &allowed);

	void triangulateSlice(const Slice &s, double z, double areaBound,
												Eigen::MatrixXd &verts, Eigen::MatrixXi &faces,
												Eigen::VectorXi &orig,
                        const std::vector<int> &allowed);

	const static double tilePadding_;

	const Slice &bottom_;
	const Slice &top_;
  Eigen::MatrixXd bbox_;
  std::vector<int> allowed_bot_, allowed_top_;
};

#endif
