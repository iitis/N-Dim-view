#pragma once

#include "dll_global.h"

#include "Mesh.h"
#include <Eigen/Geometry>

class DPVISION_DLL_API K3RoseOfWinds : public CMesh
{
public:
	K3RoseOfWinds(Eigen::MatrixXd M, double* a, CRGBA* c);
	K3RoseOfWinds(Eigen::MatrixXd M, double* a, CRGBA* c, double S);

};
