#include "K3Helpers.h"

#include "qmath.h"

#include "Mesh.h"
#include "RGBA.h"
#include "UI.h"
#include "time.h"

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
	return spectral_color(600.0 - (600.0 - 400.0) * en, a);
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

void K3ListMatrix(const wchar_t* filnam, MatrixXd XXX, const char* title) {
	int i, j;
	time_t czas;
	FILE* plik = _wfopen(filnam, L"a");
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
	// Za³ó¿my, ¿e P i P' s¹ ju¿ wype³nione danymi punktów.
	Eigen::MatrixXf P(N, 4);   // Macierz punktów przed przekszta³ceniem
	Eigen::MatrixXf P_prime(N, 4); // Macierz punktów po przekszta³ceniu

	int i, j;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < N; i++) {
			P(i, j) = Src(j, i);
			P_prime(i, j) = Dst(j, i);
		}

	}

	// Tutaj wype³nij P i P_prime odpowiednimi wartoœciami...

	// U¿yj metody najmniejszych kwadratów (SVD) do znalezienia macierzy przekszta³cenia M
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
	// std::cout << "Macierz ród³owa Src:" << std::endl << Src << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	// std::cout << "Macierz Docelowa Dst:" << std::endl << Dst << std::endl << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl << std::endl;

	// std::cout << "Macierz ród³owa P:" << std::endl << P << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	// std::cout << "Macierz Docelowa P_prime:" << std::endl << P_prime << std::endl << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl << std::endl;

	// std::cout << "Macierz przekszta³cenia M:" << std::endl << M << std::endl;
	// std::cout << "RResult:" << std::endl << RResult << std::endl;
	*/
}


void K3FillMeshToUnitCube(CMesh* ThisMesh, CRGBA Kolor) { // double l, double a) {
	CMesh* K3UnitCube = new CMesh();
	ThisMesh->addVertex(CVertex(0.0, 0.0, 0.0)); // , Kolor);
	ThisMesh->addVertex(CVertex(0.0, 0.0, 1.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(0.0, 1.0, 0.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(0.0, 1.0, 1.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(1.0, 0.0, 0.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(1.0, 0.0, 1.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(1.0, 1.0, 0.0)); //, Kolor);
	ThisMesh->addVertex(CVertex(1.0, 1.0, 1.0)); //, Kolor);
	// looking from outside, vertices should be counterclockwise:
	ThisMesh->faces().push_back(CFace(1, 5, 3));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(7, 3, 5)); // face 0

	ThisMesh->faces().push_back(CFace(4, 6, 5));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(7, 5, 6)); // face 1

	ThisMesh->faces().push_back(CFace(0, 2, 4));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(6, 4, 2)); // face 2

	ThisMesh->faces().push_back(CFace(1, 3, 0));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(2, 0, 3)); // face 3

	ThisMesh->faces().push_back(CFace(2, 3, 6));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(7, 6, 3)); // face 4

	ThisMesh->faces().push_back(CFace(0, 4, 1));
	//ThisMesh->fcolors().push_back(Kolor);
	ThisMesh->faces().push_back(CFace(5, 1, 4)); // face 5

	// UI::MESSAGEBOX::error(L"Gonna color");

	// Make front face translucent and red:
	for (int i = 0; i < 2; i++) {
		// ThisMesh->fcolors().push_back(CRGBA(0.9f, .3f, .3f, 0.3f));
		ThisMesh->fcolors().push_back(Kolor);
		// ThisMesh->fcolors().push_back(K3Color);
	};

	for (int i = 0; i < 10; i++) {
		// ThisMesh->fcolors().push_back(CRGBA(0.3f, 0.3f, 1.0f, 1.0f));
		ThisMesh->fcolors().push_back(Kolor);
		// ThisMesh->fcolors().push_back(K3Color);
	};

};


void K3LogEntry(const wchar_t* filnam, const char* mujText) {
	int i, j;
	time_t czas;
	FILE* plik = _wfopen(filnam, L"a");
	if (plik == NULL) {
		UI::MESSAGEBOX::error(L"To bardzo skomPLIKowane");
	}
	time(&czas);
	fprintf(plik, "LOG at Date=%s : %s \n", ctime(&czas), mujText);
	fclose(plik);
};



void K3ReadCSV_WithHeader(MatrixXd* X, const wchar_t* filnam) {
	char K3ParNames[20][30];

	int i;
	int j = 0;
	int nmin = 20;

	for (i = 0; i < 20; i++) {
		K3ParNames[i][0] = (char)0;
	}
	FILE* plik = _wfopen(filnam, L"r");
	if (plik) {
		int DoneRead = 8;
		j = 0;
		MatrixXd ThisSpot(1, 20);
		/* DoneRead=fscanf(plik, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s",
		K3ParNames[0], K3ParNames[1], K3ParNames[2], K3ParNames[3], K3ParNames[4],
		K3ParNames[5], K3ParNames[6], K3ParNames[7], K3ParNames[8], K3ParNames[9],
		K3ParNames[10], K3ParNames[11], K3ParNames[12], K3ParNames[13], K3ParNames[14],
		K3ParNames[15], K3ParNames[16], K3ParNames[17], K3ParNames[18], K3ParNames[19]); */
		while (DoneRead != EOF) {
			char K3buf[300];
			char* K3Yad;
			K3Yad = fgets(K3buf, 290, plik);
			if (K3Yad == K3buf) {
				DoneRead = sscanf(K3buf, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
					&ThisSpot(0, 0), &ThisSpot(0, 1), &ThisSpot(0, 2), &ThisSpot(0, 3), &ThisSpot(0, 4),
					&ThisSpot(0, 5), &ThisSpot(0, 6), &ThisSpot(0, 7), &ThisSpot(0, 8), &ThisSpot(0, 9),
					&ThisSpot(0, 10), &ThisSpot(0, 11), &ThisSpot(0, 12), &ThisSpot(0, 13), &ThisSpot(0, 14),
					&ThisSpot(0, 15), &ThisSpot(0, 16), &ThisSpot(0, 17), &ThisSpot(0, 18), &ThisSpot(0, 19));
				if (DoneRead != EOF) {
					K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", ThisSpot, "ThisSpot read");
					if ((DoneRead > 2) /* && (j < 1400) */) {  // QQ ogran.
						for (i = 0; i < DoneRead; i++) {
							(*X)(i, j) = ThisSpot(0, i);
						}
						if (DoneRead < nmin) {
							nmin = DoneRead;
						}
					}
					else
					{
						DoneRead = EOF;
					}
					j++;
				}
			}
			else
			{
				DoneRead = EOF;
			}
		}
	};
	MatrixXd A;
	A = X->block(0, 0, nmin, j);
	*X = A; // .transpose();

	fclose(plik);
};

void K3AddTotemMeshFromSpotQQQnicWeWogle(CMesh* whereto, Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook) {  // "HyperSight = DataPoint*Observer"
	CMesh* K3Korpus = new CMesh;
	Eigen::Matrix4d K3Position;
	K3Position << 1, 0, 0, K3HyperSpot[0],
		0, 2, 0, K3HyperSpot[1],
		0, 0, .6, K3HyperSpot[2],
		0, 0, 0, 1;

	/*Eigen::Matrix4d K3MakeKorpus;
	K3MakeKorpus << 1, 0, 0, 0,
	0, 2, 0, 0,
	0, 0, .2, 0,
	0, 0, 0, 1;*/

	K3FillMeshToUnitCube(K3Korpus, CRGBA(1.0f, 1.0f, 0.70f, 0.990f));  // qq kolorki

	// Scale up the torso:
	for (int i = 0; i < K3Korpus->vertices().size(); i++)
	{
		K3Korpus->vertices()[i].transformByMatrix(K3Position);
	}


	for (int i = 0; i < K3Korpus->vertices().size(); i++)
	{
		K3Korpus->vertices()[i].transformByMatrix(K3Position);
	}

	for (int i = 0; i < K3Korpus->vertices().size(); i++) {
		whereto->vertices().push_back(K3Korpus->vertices()[i]);
	};

	for (int i = 0; i < K3Korpus->faces().size(); i++) {
		whereto->faces().push_back(K3Korpus->faces()[i]);
	};

	// delete korpus;


	Eigen::Matrix4d K3MakeArm;
	K3MakeArm << 5, 0, 0, 0,
		0, .2, 0, 0,
		0, 0, .05, 0,
		0, 0, 0, 1;

	double alfa = K3HyperLook[4];

	Eigen::Matrix4d K3TwistArm;
	K3TwistArm << cos(alfa), -sin(alfa), 0, 0,
		sin(alfa), cos(alfa), 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	Eigen::Matrix4d K3ShiftArm;
	K3ShiftArm << 1, 0, 0, 2,
		0, 1, 0, 2,
		0, 0, 1, 0,
		0, 0, 0, 1;

}
