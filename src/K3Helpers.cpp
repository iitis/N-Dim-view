#include "K3Helpers.h"

#include "K3Arrow.h"
#include "Model3D.h"

#include "qmath.h"

#include "RGBA.h"
#include "UI.h"
#include "time.h"


//#define LITERAL_STRINGIFY(x) L ## x
//#define TO_WIDE_STRING(x) LITERAL_STRINGIFY(x)

QString DATA_PATH(const QString &fname)
{
	//std::wstring path(TO_WIDE_STRING(WORK_DIR));
	QString path(WORK_DIR);
	path.append(fname);
	return path;
}

//QString PNG_PATH(const QString &fname)
//{
//	QString path(PNG_DIR);
//	path.append(fname);
//	return path;
//}

//void delete_old_screenshots(const QString& pattern)
//{
//	QDir dir(PNG_DIR);
//	if (!dir.exists()) {
//		qWarning() << "Directory does not exist:" << PNG_DIR; return;
//	}
//
//	QStringList filters;
//	filters << pattern;
//	dir.setNameFilters(filters);
//	QFileInfoList fileList = dir.entryInfoList();
//	for (const QFileInfo& fileInfo : fileList) {
//		if (fileInfo.isFile()) {
//			QFile file(fileInfo.absoluteFilePath());
//			if (file.remove()) {
//				qDebug() << "Deleted:" << fileInfo.fileName();
//			}
//			else {
//				qWarning() << "Failed to delete:" << fileInfo.fileName();
//			}
//		}
//	}
//}






void K3Neg3(double A[3], double B[3]) {
	for (int i = 0; i < 3; i++) {
		B[i] = -A[i];
	}
}

CRGBA spectral_color(double l, double a) // RGB <- lambda l = < 380,780 > [nm]
{ // original:  https://stackoverflow.com/questions/22141206/how-do-i-draw-a-rainbow-in-freeglut/22149027#22149027
	double r, g, b;
	if (l < 380.0) r = 0.00;
	else if (l < 400.0) r = 0.05 - 0.05 * sin(M_PI * (l - 366.0) / 33.0);
	else if (l < 435.0) r = 0.31 * sin(M_PI * (l - 395.0) / 81.0);
	else if (l < 460.0) r = 0.31 * sin(M_PI * (l - 412.0) / 48.0);
	else if (l < 540.0) r = 0.00;
	else if (l < 590.0) r = 0.99 * sin(M_PI * (l - 540.0) / 104.0);
	else if (l < 670.0) r = 1.00 * sin(M_PI * (l - 507.0) / 182.0);
	else if (l < 730.0) r = 0.32 - 0.32 * sin(M_PI * (l - 670.0) / 128.0);
	else              r = 0.00;
	if (l < 454.0) g = 0.00;
	else if (l < 617.0) g = 0.78 * sin(M_PI * (l - 454.0) / 163.0);
	else              g = 0.00;
	if (l < 380.0) b = 0.00;
	else if (l < 400.0) b = 0.14 - 0.14 * sin(M_PI * (l - 364.0) / 35.0);
	else if (l < 445.0) b = 0.96 * sin(M_PI * (l - 395.0) / 104.0);
	else if (l < 510.0) b = 0.96 * sin(M_PI * (l - 377.0) / 133.0);
	else              b = 0.00;
	return CRGBA((float)r, (float)g, (float)b, (float)a); // (r, g, b, 255); //QQ uzmiennic alfa
}


CRGBA K3_color(double en, double a) // RGB <- energy = < 0.0, 1.0 >
{
	double l;
	if (en < 0) return CRGBA(0.0f, 0.0f, 0.0f, a);
	if (en > 1) return CRGBA(1.0f, 1.0f, 1.0f, a);
	// return spectral_color(780.0 - (780.0 - 380.0) * en, a);
	return spectral_color(600.0 - (600.0 - 400.0) * 0.7* en, a);
};


double K3ToUnity(double* X, int N) {
	// normalize N-dimensional X vector so its norm becomes unity.
	// if vector was zero, return a unity vector with all components equal.
	double s;
	double upsilon = 1e-100; //QQ
	int i;
	s = 0;
	for (i = 0; i < N; i++) {
		s += X[i] * X[i];
	};
	s = sqrt(s);
	if (s >= upsilon) {
		for (i = 0; i < N; i++) {
			X[i] /= s;
		}

	}
	else
		for (i = 0; i < N; i++) {
			X[i] = 1.0 / (double)N;
		}
	return s;
}

void K3CrossProduct(double X[3], double Y[3], double Z[3]) {
	// X, Y - given vectors;
	// result will be put in Z;
	Z[0] = X[1] * Y[2] - X[2] * Y[1];
	Z[1] = X[2] * Y[0] - X[0] * Y[2];
	Z[2] = X[0] * Y[1] - X[1] * Y[0];
};


void K3AddUnitProngs(double X[3], double Y[3], double Z[3]) {
	// create unit-length vectors perpendicular to X and to each other
	// X - given vector; results will be stored in Y and Z
	// if X={0,0,0}, results will be two random but perpendicular unit vectors

	int i, imin = 0;
	double cmin = abs(X[0]);
	for (i = 1; i < 3; i++) {
		if (abs(X[i]) < cmin) {
			imin = i;
			cmin = abs(X[i]);
		};
	}
	Y[0] = Y[1] = Y[2] = 0.0;
	Y[imin] = 1.0;
	K3CrossProduct(X, Y, Z);
	K3ToUnity(Z, 3);
	K3CrossProduct(X, Z, Y);
	K3ToUnity(Y, 3);
}

void K3ListMatrix(const QString &filnam, MatrixXd XXX, const char* title) {
	int i, j;
	time_t czas;
	FILE* plik = fopen(filnam.toStdString().c_str(), "a");
	if (plik == NULL) {
		UI::MESSAGEBOX::error(L"To bardzo skomPLIKowane");
	}
	time(&czas);
	fprintf(plik, "%s=[ %%  (%d*%d) Date=%s\n", title, XXX.rows(), XXX.cols(), ctime(&czas));
	for (int j = 0; j < XXX.rows(); j++) {
		for (i = 0; i < XXX.cols(); i++)
		{
			if (i > 0) {
				fprintf(plik, ", ");
			};
			fprintf(plik, "%15lg", XXX(j, i));
		}
		fprintf(plik, "\n");
		//		fprintf(plik, "%15lg,  %15lg,  %15lg\n", XXX(j, 0), XXX(j, 1), XXX(j, 2));
	}
	fprintf(plik, "]\n");
	fclose(plik);
};

// Eigen::Matrix4d K3FindTransform(Matrix4Xd Src(4, N), Matrix4Xd Dst(4, N), int N);


void K3FindTransform(Eigen::Matrix4d Src, Eigen::Matrix4d Dst, Eigen::Matrix4d* RResult, int N) {
	// Z maila od Darka - dziekuje! QQ
	// Załóżmy, że P i P' są już wypełnione danymi punktów.
	Eigen::MatrixXf P(N, 4);   // Macierz punktów przed przekształceniem
	Eigen::MatrixXf P_prime(N, 4); // Macierz punktów po przekształceniu

	int i, j;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < N; i++) {
			P(i, j) = Src(j, i);
			P_prime(i, j) = Dst(j, i);
		}

	}

	// Tutaj wypełnij P i P_prime odpowiednimi wartościami...

	// Użyj metody najmniejszych kwadratów (SVD) do znalezienia macierzy przekształcenia M
	Eigen::MatrixXf M = P.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(P_prime);


	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
			(*RResult)(i, j) = M(j, i);
			(*RResult)(i, j) += (i == j ? 0.1 : 0.0);

		}
		// (*RResult)(j, j) +=0.1;
	}

	// std::cout << "X";
	/*
	// std::cout << "Macierz Źródłowa Src:" << std::endl << Src << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	// std::cout << "Macierz Docelowa Dst:" << std::endl << Dst << std::endl << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl << std::endl;

	// std::cout << "Macierz Źródłowa P:" << std::endl << P << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	// std::cout << "Macierz Docelowa P_prime:" << std::endl << P_prime << std::endl << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl << std::endl;

	// std::cout << "Macierz przekształcenia M:" << std::endl << M << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	*/
}



void K3FillMat(MatrixXd& X, double a) {
	int j, i;
	for (j = 0; j < X.rows(); j++) {
		for (i = 0; i < X.cols(); i++) {
			X(j, i) = a;
		}
	}
}


void K3_4x4viewN(MatrixXd* V, int k, double alfa) {

	/*
0, 1, 2, 3,
0, 2, 1, 3,
0, 3, 1, 2,
1, 3, 0, 2,
1, 2, 0, 3,
3, 2, 0, 1*/

	int k3tasOLD[6][4] =
	{
0, 1, 2, 3,  // 1 2 swapped
0, 2, -1, 3,  // 1 3
0, 3, -1, -2,  // 0 2DR
-1, -2, -0, -3,  // 0 3
-3, -2, -0, 1 };
	int k3tas[6][4] =
	{
		1, 2, 3, 4,  // 1 2 swapped
			1, 3, -2, 4,  // 1 3
			1, 4, -2, -3,  // 0 2
			-2, 4, -1, -3,  // 1 3
			-2, -3, -1, -4,  // 0 3
			-4, -3, -1, 2 };
	int k3_TurnPlane[6][2] = {
		1,2,
		1,3,
		0,2,
		1,3,
		0,3,
		0,2 };
	K3FillMat(*V, 0.0);
	for (int i = 0; i < 4; i++) {
		if (k3tas[k][i] > 0) {
			(*V)(i, abs(k3tas[k][i]) - 1) = 1.0;
		}
		else {
			(*V)(i, abs(k3tas[k][i]) - 1) = -1.0;
		}
	}
	// k = k % 3;
	int ia = k3_TurnPlane[k][0];
	int ib = k3_TurnPlane[k][1];

	double sa = sin(alfa);
	double ca = cos(alfa);
	double a4[4] = { 0, 0, 0, 0 };
	double b4[4] = { 0, 0, 0, 0 };
	for (int i = 0; i < 4; i++) {
		a4[i] += ca * (*V)(ia, i) + sa * (*V)(ib, i);
		b4[i] += ca * (*V)(ib, i) - sa * (*V)(ia, i);
	}
	for (int i = 0; i < 4; i++) {
		(*V)(ia, i) = a4[i];
		(*V)(ib, i) = b4[i];
	};
}



void K3ArrowsArc(double Center[3], double A[3], double B[3], CModel3D* K3MyModel, double R, int n, CRGBA* colour) {
	// make arc from Center+A to Center+B, made of 12 short arrows
	int i, i88, i99;
	double a[3], b[3];
	for (int i = 0; i < n + 0; i++) {
		for (int i88 = 0; i88 < 3; i88++) {
			double Rup = 1.0 + 1.0 * ((double)(2 * n - 2 * i)) * ((double)(2 * i)) / ((double)(n * n));
			a[i88] = A[i88] * Rup * ((double)(n - i) / (double)n) + B[i88] * Rup * ((double)(i) / (double)n);
			Rup = 1.0 + 1.0 * ((double)(2 * n - 2 * (i + 1))) * ((double)(2 * (i + 1))) / ((double)(n * n));
			b[i88] = A[i88] * Rup * ((double)(n - i - 1) / (double)n) + B[i88] * Rup * ((double)(i + 1) / (double)n);
			double da;
			da = 0.3 * (b[i88] - a[i88]);
			b[i88] += da / 3.0;
			a[i88] -= da / 3.0;
			a[i88] += Center[i88];
			b[i88] += Center[i88];

		};
		K3Arrow* Arro1 = new K3Arrow(a, b, R, colour);
		K3MyModel->addChild(Arro1);
	}
};


