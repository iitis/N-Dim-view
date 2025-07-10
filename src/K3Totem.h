#pragma once

#include "dll_global.h"

#include "Mesh.h"
#include "K3ChernoffFace.h"

#include<Eigen/Geometry>
#include <optional>

class DPVISION_DLL_API K3Totem : public CMesh { // CModel3D {
	// extra features will be added here

public:
	// K3Totem(std::vector<double> K3HyperSpot);  // "K3HyperSpot = DataPoint*Observer"
	K3Totem(Eigen::VectorXd P, Eigen::VectorXd V);  // "K3HyperSpot = DataPoint*Observer"

	void K3FillMeshToUnitCube(CMesh* ThisMesh, CRGBA Kolor);

	std::shared_ptr<K3ChernoffFace> make_face(std::vector<std::optional<double>> values, CMesh* korpus);
	
	// AdjustData(std::vector<double> K3HyperSpot);  Eigen::VectorXd
};

