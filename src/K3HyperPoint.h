#pragma once

#include <Eigen/Geometry>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "dll_global.h"

class DPVISION_DLL_API K3HyperPoint {
public:
	std::vector<double> coords;
	K3HyperPoint(int n = 3);

	void setN(int n);
	//void K3ShowTotem();
};


