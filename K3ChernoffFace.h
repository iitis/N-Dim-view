#pragma once

#include "dll_global.h"

#include <Image.h>

class DPVISION_DLL_API K3ChernoffFace : public CImage {
public:
	K3ChernoffFace(int cx, int cy, double k3params[10]);
};
