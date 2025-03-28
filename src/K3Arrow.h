#pragma once

#include "dll_global.h"

#include "Mesh.h"

class CRGBA;

class DPVISION_DLL_API K3Arrow : public CMesh {
public:
	K3Arrow(double A[3], double B[3], double R, CRGBA* colour);
};

