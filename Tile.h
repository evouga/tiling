#ifndef TILE_H
#define TILE_h

#include <Eigen/Core>

class Slice;

class Tile
{
public:
	Tile(const Slice &bottom, const Slice &top);
	~Tile();
	void triangulateSlices(double areaBound, 
		Eigen::MatrixXd &botverts, Eigen::MatrixXi &botfaces, 
		Eigen::MatrixXd &topverts, Eigen::MatrixXi &topfaces);

	

private:
	void computeCubeTransformation();
	void triangulateSlice(const Slice &s, double z, double areaBound, Eigen::MatrixXd &verts, Eigen::MatrixXi &faces);

	const static double tilePadding_;

	const Slice &bottom_;
	const Slice &top_;	
	Eigen::Vector2d translate_;
	Eigen::Matrix2d scale_;
};

#endif
