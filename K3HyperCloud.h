#pragma once

#include "dll_global.h"

#include <Eigen/Geometry>

using namespace Eigen;

#include "K3HyperPoint.h"

class DPVISION_DLL_API K3HyperCloud { // possibly : public Eigen::MatrixXd {
public:
	std::vector<K3HyperPoint> K3RawData;    // QQ <K3HyperPoint>
	// std::vector <std::vector <double, 3>, 15> Observer;   // dimensions temporary
	MatrixXd Observer = MatrixXd::Constant(15, 15, 0.0);
	// std::vector <double[3][3]> Observer;   // dimensions temporary
	int n, m; // n rows, m columns;

	// K3RawData[m][n]  dostęp do n-tej wsp. m-tego punktu

	K3HyperCloud(int new_n = 4, int new_m = 4);
};


