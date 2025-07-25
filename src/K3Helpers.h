#pragma once

#include "dll_global.h"

#include<Eigen/Geometry>
using namespace Eigen;

class CRGBA;
class CMesh;
class CModel3D;

class QString;

#define WORK_DIR "c:/K3/Wielowymiar/"

QString DPVISION_DLL_API DATA_PATH(const QString &fname);

void DPVISION_DLL_API K3Neg3(double A[3], double B[3]);
CRGBA DPVISION_DLL_API K3_color(double en, double a);
double DPVISION_DLL_API K3ToUnity(double* X, int N);
void DPVISION_DLL_API K3AddUnitProngs(double X[3], double Y[3], double Z[3]);
void DPVISION_DLL_API K3ListMatrix(const QString &filnam, MatrixXd XXX, const char* title);
void DPVISION_DLL_API K3FindTransform(Eigen::Matrix4d Src, Eigen::Matrix4d Dst, Eigen::Matrix4d* RResult, int N);

void DPVISION_DLL_API K3FillMat(MatrixXd& X, double a);
void DPVISION_DLL_API K3_4x4viewN(MatrixXd* V, int k, double alfa);

void DPVISION_DLL_API get_observer_matrix(Eigen::MatrixXd& dst, int k, double alfa);
Eigen::Array<bool, Eigen::Dynamic, 1> DPVISION_DLL_API create_slab_mask(Eigen::MatrixXd& V, Eigen::MatrixXd& X_spatial, double slab_threshold);
Eigen::MatrixXd DPVISION_DLL_API use_mask(Eigen::MatrixXd& X_view, Eigen::Array<bool, Eigen::Dynamic, 1>& mask);

void DPVISION_DLL_API K3ArrowsArc(double Center[3], double A[3], double B[3], std::shared_ptr<CModel3D> K3MyModel, double R, int n, CRGBA* colour);
