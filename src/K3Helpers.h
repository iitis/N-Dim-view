#pragma once

#include "dll_global.h"

#include<Eigen/Geometry>
using namespace Eigen;

class CRGBA;
class CMesh;
class CModel3D;


#define WORK_DIR "c:/K3/Wielowymiar/"
#define PNG_DIR "d:/K3/Wielowymiar/"

QString DPVISION_DLL_API DATA_PATH(const QString &fname);
QString DPVISION_DLL_API PNG_PATH(const QString &fname);
void DPVISION_DLL_API delete_old_screenshots(const QString& pattern);

void DPVISION_DLL_API K3Neg3(double A[3], double B[3]);
CRGBA DPVISION_DLL_API K3_color(double en, double a);
double DPVISION_DLL_API K3ToUnity(double* X, int N);
void DPVISION_DLL_API K3AddUnitProngs(double X[3], double Y[3], double Z[3]);
void DPVISION_DLL_API K3ListMatrix(const QString &filnam, MatrixXd XXX, const char* title);
void DPVISION_DLL_API K3FindTransform(Eigen::Matrix4d Src, Eigen::Matrix4d Dst, Eigen::Matrix4d* RResult, int N);

void DPVISION_DLL_API K3FillMat(MatrixXd& X, double a);
void DPVISION_DLL_API K3_4x4viewN(MatrixXd* V, int k, double alfa);

void DPVISION_DLL_API K3ArrowsArc(double Center[3], double A[3], double B[3], CModel3D* K3MyModel, double R, int n, CRGBA* colour);