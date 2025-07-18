#pragma once

#include "dll_global.h"

#include "Mesh.h"
#include "K3ChernoffFace.h"

#include<Eigen/Geometry>
#include <optional>
#include <memory>

class DPVISION_DLL_API K3Totem : public CMesh { // CModel3D {
	// extra features will be added here

public:
	// K3Totem(std::vector<double> K3HyperSpot);  // "K3HyperSpot = DataPoint*Observer"
	//K3Totem(Eigen::VectorXd P, Eigen::VectorXd V);  // "K3HyperSpot = DataPoint*Observer"

	static std::shared_ptr<K3Totem> create(Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook);

	void K3FillMeshToUnitCube(std::shared_ptr<CMesh> ThisMesh, CRGBA Kolor);

	std::shared_ptr<K3ChernoffFace> make_face(std::vector<std::optional<double>> values, std::shared_ptr<CMesh> korpus);
	
	// AdjustData(std::vector<double> K3HyperSpot);  Eigen::VectorXd
};

