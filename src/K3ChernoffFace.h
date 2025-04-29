#pragma once

#include "dll_global.h"

#include <Image.h>
#include <optional>

class DPVISION_DLL_API K3ChernoffFace : public CImage {
public:
	K3ChernoffFace(int cx, int cy, std::vector<std::optional<double>> k3params);
};
