#pragma once

#include "dll_global.h"

#include<Eigen/Geometry>
using namespace Eigen;

class CRGBA;
class CMesh;

CRGBA DPVISION_DLL_API K3_color(double en, double a);
double DPVISION_DLL_API K3ToUnity(double* X, int N);
void DPVISION_DLL_API K3AddUnitProngs(double X[3], double Y[3], double Z[3]);
void DPVISION_DLL_API K3ListMatrix(const wchar_t* filnam, MatrixXd XXX, const char* title);
void DPVISION_DLL_API K3FillMeshToUnitCube(CMesh* ThisMesh, CRGBA Kolor);
void DPVISION_DLL_API K3FindTransform(Eigen::Matrix4d Src, Eigen::Matrix4d Dst, Eigen::Matrix4d* RResult, int N);

