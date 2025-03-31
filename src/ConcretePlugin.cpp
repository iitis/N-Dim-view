/*
MISJA:
Stworzyć nawigację nD

- Założenia: nie korzystam z ruszania myszką, wszystko programowo.

- czy można reagować na klawisze?

- modyfikacja macierzy położenia obserwatora:

  -  sterowanie 3D w aktualnej przestrzeni obserwacji;
  - jak się obracać ku innym wymiarom?
	 - do kolumny "ku" dodawać porcjami,
		   w wierszach x? y? z?
  - Jak wyświetlać n-wymiarową różę wiatrów?
	- strzałki z wersorów tych wymiarów, które mają udział nie mniejszy niż...

  -  skrypty w plikach tekstowych;

  CEL NAJBLIŻSZY:
  Trajektoria obserwatora, założenia na pozycjonowanie w nD.

  - 3 widzialne wymiary ustalać w kolejności zgodnej z numeracją w hiperprzestrzeni
  - A jeżeli są mieszane? W kolejności malejącej średniej ważonej indeksów w nD.

  Zakładam - przynajmniej na razie - że wymiary odzwierciedlane wyglądem
	  (a nie lokalizają) totemu będą "pojedyncze", tzn każda cecha wyglądu
	  oddaje tylko jeden parametr źródłowy. To ważne ograniczenie na macierz Observer.

  Ruchy 3 widzialnych wymiarów w przestrzeni nD:
  PILOTEM nazywam źródło komend dla obserwatora, może nim być:
	  użytkownik przy klawiaturze,
	  skrypt wczytywany z pliku albo z tablicy stałych
  Komendą może być:
	  identyfikator (numer? nazwa?) któregoś wymiaru z nD;
		  taką komendę wykonuje sie przez zwiększenie udziału wskazanego wymiaru...
			- ale w czym i kosztem czego? Pasowałoby w tym z Trzech,
			   który już _tego_ wymiaru ma najwięcej. Jeśli żaden nie ma, to w tym,
			   który jest mu najbliżej średnią ważoną indeksów.
			   I kosztem dotychczasowych składowych tegoż wymiaru.
			   Tylko: solidarnie czy najchudszy traci?
	  identyfikator (numer? nazwa?) któregoś wymiaru z 3D;
		  taką komendę wykonuje sie przez obrót widzialnego przekroju
			 wokół wskazanego wymiaru o 45 stopni


   PROGRAM PRACY:
   + napisać funkcję, wyznaczającą średni (ważenie) indeks w nD danego wymiaru 3D:
   K3AverIndex(Observer, iSmallDim);
   + Funkcja regulująca udział danego (indeksem) wymiaru nD w danym (indeksem) w. 3D
   + K3AdjustShare



*/

#include "ConcretePlugin.h"

#include "AP.h"

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues> 

#include "K3Helpers.h"

#include "K3Totem.h"
#include "K3Arrow.h"
#include "K3ChernoffFace.h"
#include "K3RoseOfWinds.h"

#include "CsvReader.h"

using namespace Eigen;

// dodatek LL:
void K3LogEntry(const wchar_t* filnam, const char* myText);

#define K3GRAND_SCALE 50.0
// było 0.5
// QQ Scale!
// #include<Eigen/Geometry>
// #include <Eigen/SVD> // QQ czy dobrze?
// number of human-identifiable parameters to be displayed hyperspacially:
#define K3hs 2  // density, pH
// number of human-identifiable parameters to be displayed by avatar appearance:
#define K3hv 2  // citric, suger
// number of NOT human-identifiable parameters to be displayed hyperspacially:
#define K3as 3
// reduces to:
#define K3asr 2
// number of NOT human-identifiable parameters to be displayed by avatar appearance:
#define K3av 6
// reduces to:
#define K3avr 4

// ============================================================
// (DP): FUNKCJE DODANE PRZEZE MNIE W CELU UJEDNOLICENIA
// I UPROSZCZENIA POTENCJALNEJ ZMIANY ŚCIEŻEK ZAPISU
// Proponuje w kodzie zamiast pełnej ścieżki używać: DATA_PATH(<nazwa_pliku>)
// Możesz mieć oddzielny katalog dla danych i dla wynikowych screenshotów

#define WORK_DIR "c:/K3/Wielowymiar/"
#define PNG_DIR "d:/K3/Wielowymiar/"

#define LITERAL_STRINGIFY(x) L ## x
#define TO_WIDE_STRING(x) LITERAL_STRINGIFY(x)

std::wstring DATA_PATH(std::wstring fname)
{
	std::wstring path(TO_WIDE_STRING(WORK_DIR));
	path.append(fname);
	return path;
}

QString PNG_PATH(QString fname)
{
	QString path(PNG_DIR);
	path.append(fname);
	return path;
}

void delete_old_screenshots(const QString& pattern)
{
	QDir dir(PNG_DIR);
	if (!dir.exists()) {
		qWarning() << "Directory does not exist:" << PNG_DIR; return;
	}

	QStringList filters;
	filters << pattern;
	dir.setNameFilters(filters);
	QFileInfoList fileList = dir.entryInfoList();
	for (const QFileInfo& fileInfo : fileList) {
		if (fileInfo.isFile()) {
			QFile file(fileInfo.absoluteFilePath());
			if (file.remove()) {
				qDebug() << "Deleted:" << fileInfo.fileName();
			}
			else {
				qWarning() << "Failed to delete:" << fileInfo.fileName();
			}
		}
	}
}

// ============================================================

//#include "Image.h"


void K3FillMat(MatrixXd& X, double a) {
	int j, i;
	for (j = 0; j < X.rows(); j++) {
		for (i = 0; i < X.cols(); i++) {
			X(j, i) = a;
		}
	}
}



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


void OldK3_4x4viewN(MatrixXd* V, int k, double alfa) {

	int k3tas[13][4] =
	{ 0, 1,	2,	3,  // 01
	1,	0,	2,	3,  // 12
	1,	2,	0,	3,  // 30
	3,	2,	0,	1,  // 01
	2,	3,	0,	1,  // 12
	2,	0,	3,	1,  // 30
	1,	0,	3,	2,  // 01
	0,	1,	3,	2,  // 12
	0,	3,	1,	2,  // 30
	2,	3,	1,	0,  // 01
	3,	2,	1,	0,  // 12
	3,	1,	2,	0,  // 30
	0,	1,	2,	3 };// 01
	int k3_TurnPlane[3][2] = {
		0,1,
		1,2,
		3,0 };
	K3FillMat(*V, 0.0);
	for (int i = 0; i < 4; i++) {
		(*V)(i, k3tas[k][i]) = 1.0;
	}
	k = k % 3;
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


double K3AverIndex(MatrixXd Observer, int iSmallDim) {
	double SumNumerator, SumDenominator, A;
	SumNumerator = SumDenominator = 0.0;
	int i, j;
	for (i = 0; i < Observer.cols(); i++) {
		SumNumerator += Observer(iSmallDim, i);
		SumDenominator += (double)i * Observer(iSmallDim, i);
	}
	if (abs(SumDenominator) > 0.001) {
		A = SumNumerator / SumDenominator;
	}
	else {
		A = 0.0;
	};
	return A;
}

void K3AdjustShare(MatrixXd Observer, int iBigDim, int iSmallDim, double delta) {
	// QQ Uzupelnic komentarz - co robi ta funkcja??
	int i;
	double SumOfAbsolutes, Factor;
	SumOfAbsolutes = 0.0;
	for (i = 0; i < Observer.cols(); i++) { // QQ może warto pominąć iBigDim?
		SumOfAbsolutes += abs(Observer(iSmallDim, i));
	};
	if (SumOfAbsolutes > 0.001) {
		Factor = (SumOfAbsolutes - abs(delta)) / SumOfAbsolutes;
		if (Factor < 0.0) {
			Factor = 0.0;
		}
		for (i = 0; i < Observer.cols(); i++) {  // QQ i tutaj też?
			Observer(iSmallDim, i) *= Factor;
		};

	};
	Observer(iSmallDim, iBigDim) += delta;
};




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
		K3Arrow* Arro1 = new	K3Arrow(a, b, R, colour);
		K3MyModel->addChild(Arro1);

	}
};


// Now to the job:

#define K3Dimensionality 15

// :koniec dodatku LL


ConcretePlugin::ConcretePlugin(void)
{
	using namespace Eigen;

	//UI::MENU::addUserAction( "Test menu", myAction );
	//UI::MESSAGEBOX::error( L"I'm LOADED" );

	m_picking = false;

	current_data_matrix = Eigen::MatrixXd();
	current_assignment = QVector<ColumnAssignment>();

	//MatrixXd X77(200, 6);
	//K3ReadCSV_WithHeader(&X77, L"C:/K3/Wielowymiar/Tabulka45.txt");
	//K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", X77, "X77 iii");

	// K3 loading CSV file:
	char K3ParNames[12][30];
	int j;
	// char K3Bufor[80];
	Eigen::Matrix3d m3;
	m3(1, 2) = 7112;


	// K3 trying matrices:
	// m3(0, 0) =   // that was in the original, let's try to be compact:
	m3(1, 2) = 713;
	// m3 << 0, 1, 2, 10, 11, 12, 21, 22, 23;
	// :end of K3 trying matrices

}

ConcretePlugin::~ConcretePlugin(void)
{
}


void ConcretePlugin::onLoad()
{
	UI::PLUGINPANEL::create(m_ID, L"Wielowymiar");

	UI::PLUGINPANEL::addButton(m_ID, L"csvReaderTest", L"Wczytaj plik CSV", 0, 0);

	UI::PLUGINPANEL::addButton(m_ID, L"K3Display", L"K3Display", 6, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3CzteroPajak", L"Zrób K3CzteroPajak", 7, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3Krata", L"Zrób Kratę", 8, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3Dance", L"K3Dance", 9, 0);
}


//#include "AnnotationPlane.h"
//#include "AnnotationPoint.h"
//#include <vector>
//#include <Eigen/src/Eigenvalues/EigenSolver.h>



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
// {
//	CImage* im = new CImage(cx, cy, CImage::Format::Format_ARGB32);


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

};


// Principal Component Analysis:
MatrixXd ConcretePlugin::K3_Get_PCA_Funnel(MatrixXd X, int nd) {
	// Given an n*m matrix X, and an integer 0<nd<n,
	// create a result nd*m matrix R such that
	// muliplying R*X will yield the data in reduced dim.
	int n, m, i, j;
	n = X.rows();
	m = X.cols();
	MatrixXd R(nd, n);
	// QQ do something!

	// Compute row averages RA:
	MatrixXd RA(n, 1);
	MatrixXd Adder(m, 1);
	K3FillMat(Adder, 1.0 / ((double)m));
	RA = X * Adder;
	// get B = centered copy of X:
	MatrixXd B = X;  // see Wikipedia on PCA
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			B(j, i) -= RA(j, 0);
		};
	};
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", B, "B");

	BDCSVD<MatrixXd> svd(B, ComputeFullU | ComputeFullV);

	MatrixXd  U = svd.matrixU();  // used to be auto
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", U, "U");
	MatrixXd V = svd.matrixV();
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", V, "V");
	MatrixXd sigma = svd.singularValues().asDiagonal().toDenseMatrix();
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", sigma, "sigma");

	return U.block(0, 0, nd, n);
};

// Principal Component Analysis:
MatrixXd /* ConcretePlugin:: */ oldK3_Get_PCA_Funnel(MatrixXd X, int nd) {
	// Given an n*m matrix X, and an integer 0<nd<n,
	// create a result nd*m matrix R such that
	// muliplying R*X will yield the data in reduced dim.

	int n, m, i, j;
	n = X.rows();
	m = X.cols();
	MatrixXd R(nd, n);
	// QQ do something!

	// Compute row averages RA:
	MatrixXd RA(n, 1);
	MatrixXd Adder(m, 1);
	K3FillMat(Adder, 1.0 / ((double)m));
	RA = X * Adder;
	MatrixXd B = X;  // see Wikipedia on PCA
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			B(j, i) -= RA(j, 0);
		}

	}

	// Now compute covariance:
	MatrixXd Covar = B * B.transpose();

	// Now scale B according to Covar:
	MatrixXd BScal = B;
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			BScal(j, i) /= sqrt(Covar(j, j) + 0.00001);
		}

	}

	// Now compute correlation:
	MatrixXd Corr = BScal * BScal.transpose();

	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", Corr, "Do rozkladu Corr");


	// https://blog.demofox.org/2022/07/12/calculating-svd-and-pca-in-c/
		// https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca


		// Matrix2d CTC = C.transpose() * C;

		//EigenSolver<Matrix2d> es(Corr);
		//EigenSolver<Matrix2d>::EigenvalueType eigenValues = es.eigenvalues();
		//EigenSolver<Matrix2d>::EigenvectorsType eigenVectors = es.eigenvectors();

	BDCSVD<Matrix2d> svd(Corr, ComputeFullU | ComputeFullV);

	MatrixXd  U = svd.matrixU();  // used to be auto
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", U, "U");
	MatrixXd V = svd.matrixV();
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", V, "V");
	MatrixXd sigma = svd.singularValues().asDiagonal().toDenseMatrix();
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", sigma, "sigma");

	// vhttps://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html#a9e03d09fd7cfc0120847fcb63aa353f3

	EigenSolver<MatrixXd> K3es;
	K3es.compute(Corr, /* computeEigenvectors = */ true);
	auto K3EV = K3es.eigenvectors();
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3EV, "K3EV");;

	for (j = 0; j < nd; j++) {
		for (i = 1; i < n; i++) {
			R(j, i) = K3EV(i, j).real();
		}
	}
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", R, "Finishing R");

	return R;
};

void  ConcretePlugin::K3Display(CModel3D* o, float d)
{
	CMesh* dstMesh = (CMesh*)(o->getChild());

	dstMesh->calcVN();

	for (int i = 0; i < dstMesh->vnormals().size(); i++)
	{
		CVertex* v1 = &dstMesh->vertices()[i];
		CVector3f* vn1 = &dstMesh->vnormals()[i];

		v1->Set(v1->X() + d * vn1->X(), v1->Y() + d * vn1->Y(), v1->Z() + d * vn1->Z());
		AP::processEvents();
	}
}

// RotationMat = ProjectionMat.block<3, 3>(0, 0);
// translationVec = ProjectionMat.block<3, 1>(0, 3);


/*   COPY:
		// number of human-identifiable parameters to be displayed hyperspacially:
#define K3hs 2  // density, pH
// number of human-identifiable parameters to be displayed by avatar appearance:
#define K3hv 2  // citric, suger
// number of NOT human-identifiable parameters to be displayed hyperspacially:
#define K3as 3
// reduces to:
#define K3asr 2
// number of NOT human-identifiable parameters to be displayed by avatar appearance:
#define K3av 6
// reduces to:
#define K3avr 4 */

void DPVISION_DLL_API K3AddMyCloud(CModel3D* K3MyModel, MatrixXd K3ObsCloud, MatrixXd K3ViewMat, double K3Toler) {
	int k = K3ObsCloud.cols();
	if (k > 200) {
		k = 200;
	};
	int m = K3ObsCloud.rows();
	MatrixXd K3LocalSpots(5, k);  // QQ: FUTURE: 5 = X Y Z d T (3D + distance from current 3D space + homogeneous extra dim)
			// RIGHT??: 5 = X Y Z T d (4D + distance from current 3D space)
	for (int i = 0; i < k; i += 10) {  // DARKU TUTAJ DALEM 10, ZEBY WYSWIETLAC MNIEJ ELEMENTOW
		for (int j = 0; j < 5; j++) {
			// std::cout << "shame on i, j, k:" << i << " " << j << " " << k << std::endl;
			if (j < m) {  // missing dimensions? Fill with ones.
				K3LocalSpots(j, i) = K3ObsCloud(j, i);
			}
			else {
				K3LocalSpots(j, i) = 1.0;
			}
		}
		K3LocalSpots(4, i) = 1.0;  // homogeneous extra dim  // K3UNBLOCKED 2024.08.22
					   // QQ7  was "k", replaced with "i"!
	}
	// K3ListMatrix(L"C:/K3/Wielowymiar/M1142", K3ViewMat, "K3ViewMat*");
	// K3ListMatrix(L"C:/K3/Wielowymiar/M1142", K3LocalSpots, "*K3LocalSpots");

	K3LocalSpots = K3ViewMat * K3LocalSpots; // K3ObsCloud.block(0, 0, K3hs + K3asr + 1, k); //K3hs + K3asr, k>

	// K3ListMatrix(L"C:/K3/Wielowymiar/M1142", K3LocalSpots, "=K3LocalSpots");
	// QQSize  ViewMat applied here
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3ViewMat, "K3ViewMat 1246");
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3LocalSpots, "K3LocalSpots 1246");

	for (int i = 0; i < k; i++) {
		// if ((abs(K3LocalSpots(3, i)) < K3Toler) && (abs(K3LocalSpots(4, i)) != 0))
		{
			K3Totem* K3TenTotem;
			i = i; //qq
			Eigen::VectorXd P(4);
			Eigen::VectorXd V(m - 4);
			// QQSize P is the position vector!
			double askala = 0.1 / K3LocalSpots(4, i);
			P << K3LocalSpots(0, i) * askala, K3LocalSpots(1, i)* askala, K3LocalSpots(2, i)* askala, K3LocalSpots(3, i)* askala; // Totem
			// P << 15.0 * i, 15.0 * i + 0.1, 15.0 * i + 0.2, 4.04; // Totem
			for (int j88 = 0; j88 < m - 4; j88++) {
				V(j88) = K3ObsCloud(4 + j88, i);
			}
			// K3ListMatrix(const wchar_t* filnam, MatrixXd XXX, const char* title)
			// K3ListMatrix(L"C:/K3/Wielowymiar/M1120", K3ObsCloud, "M1120");
			double K3P77[4];
			K3P77[0] = P(0);
			K3P77[1] = P(1);
			K3P77[2] = P(2);
			K3P77[3] = P(3);
			// QQSize Creating the totem!
			// K3LogEntry(L"C:/K3/Wielowymiar/MujZrzut.txt", "TotemDone");
			K3TenTotem = new K3Totem(P, V); // K3ObsCloud.block(5, i, K3hv + K3avr, 1));
			// K3ListMatrix(L"C:/K3/Wielowymiar/M1332.txt", K3ObsCloud, "TotemDone");
			K3MyModel->addChild(K3TenTotem); // QQ blok Wlkanoc
			// K3ListMatrix(L"C:/K3/Wielowymiar/M1120", K3ObsCloud, "TotemAdded");

			if (i % 500 == 0) {
				UI::updateAllViews();
				if (i % 1000 == 0) {
					std::cout << i << "... ";
				}
			}
		}
	}

	// UI::updateAllViews();
	std::cout << k << std::endl;

	//	K3ListMatrix(DATA_PATH(L"M1120").c_str(), K3ObsCloud, "CloudAdded");
}

void K3Neg3(double A[3], double B[3]) {
	for (int i = 0; i < 3; i++) {
		B[i] = -A[i];
	}
}

int ConcretePlugin::K3FormProjectionMatrix(Eigen::MatrixXd* RawData) {
	// Use list of parameter choices and data matrix to create a projection matrix
	// that projects the raw dimensions into the observables according to user wishes
	int K3I_ListAnon[100];
	int i88, j88, j_raw, j;
	int K3I_taken_s[4] = { 0,0,0,0 };
	int K3I_taken_v[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int K3I_Count_a = 0, K3I_Count_s = 0, K3I_Count_v = 0;
	if (RawData->rows() * RawData->cols() > 0)
	{
		int K3_IDim = RawData->rows();
		int K3_ICard = RawData->cols();
		if (RawData != NULL) {
		};
		K3_IObs = new Eigen::MatrixXd(4 + 10, K3_IDim); // 4 spatial + 10 face params.
		K3FillMat(*K3_IObs, 0.0);

		for (int j = 0; j < K3_IDim; j++) {
			int j_raw = current_assignment[j].featureIndex;
			int j_boiled = 0;
			switch (current_assignment[j].groupIndex) {
			case -1:  // ignore
				break;
			case 0:  // unnamed
				K3I_ListAnon[K3I_Count_a++] = j_raw;
				break;
			case 1:  // spacial
				if (current_assignment[j].label_id.has_value()) {
					j_boiled = current_assignment[j].label_id.value();
					K3I_taken_s[j_boiled] = 2;
					K3I_Count_s++;
				}
				(*K3_IObs)(j_boiled, current_assignment[j].featureIndex) = 1;
				break;
			case 2:  // visual
				if (current_assignment[j].label_id.has_value()) {
					j_boiled = current_assignment[j].label_id.value() + 4;
					K3I_Count_v++;
				};
				(*K3_IObs)(j_boiled, current_assignment[j].featureIndex) = 1;
				break;

			}
		};
		// make funnel. Firt, get the anonymous together:
		MatrixXd X_anon(K3I_Count_a, K3_ICard);
		for (j88 = 0; j88 < K3I_Count_a; j88++) {
			j_raw = K3I_ListAnon[j88];
			for (i88 = 0; i88 < K3_ICard; i88++) {
				X_anon(j88, i88) = (*RawData)(j_raw, i88);
			};
		};
		if (K3I_Count_s < 4) {
			int j_boiled = 0;
			MatrixXd X88 = K3_Get_PCA_Funnel(X_anon, 4 - K3I_Count_s); // QQ K3I_Count_a);
			K3ListMatrix(DATA_PATH(L"MojeLeje.txt").c_str(), X88, "Lejek");
			for (j88 = 0; j88 < 4 - K3I_Count_s; j88++) {
				while (K3I_taken_s[j_boiled] > 0) { // skip the takien rows
					j_boiled++;
				}; // qq uwaga na przekroczenie
				if (j_boiled < 4) {
					K3I_taken_s[j_boiled] = 3;
					for (i88 = 0; i88 < K3I_Count_a; i88++) {
						(*K3_IObs)(j_boiled, i88) = X88(j88, i88);
					};
				};

			};

		};

		// K3_IObs READY!

		// Now standardize the sigmas:

		K3BoiledData = *K3_IObs * *RawData;
		K3ListMatrix(DATA_PATH(L"MujStat.txt").c_str(), *RawData, "RawData");
		K3ListMatrix(DATA_PATH(L"MujStat.txt").c_str(), *K3_IObs, "K3_IObs");
		K3ListMatrix(DATA_PATH(L"MujStat.txt").c_str(), K3BoiledData, "K3BoiledData");
		for (j = 0; j < K3BoiledData.rows(); j++) {
			double S1 = 0.0, S2 = 0.0, M, D;
			int N = K3BoiledData.cols()-2, i;
			for (i = 1; i < N; i++) {
				S1 += K3BoiledData(j, i);
				S2 += (K3BoiledData(j, i) * K3BoiledData(j, i));
			};
			M = S1 / (double)N;
			D = sqrt(S2 / (double)N - M * M);
			if (D < 0.0001) {  // avoid division by zero
				D = 0.0001;
			};
			for (i = 1; i < N; i++) {
				// **K3_IObs.
				K3BoiledData(j, i) = (K3BoiledData(j, i) -M)/ D ;
				if (j < 4) {
					K3BoiledData(j, i) = K3BoiledData(j, i) + 0.5;
				};
			};
		};

	};
	K3ListMatrix(DATA_PATH(L"MujStat.txt").c_str(), *K3_IObs, "NEW_K3_IObs");

	return(1);
};
void ConcretePlugin::onButton(std::wstring name)
{
	//Eigen::MatrixXd *K3_IObs;
	if (0 == name.compare(L"csvReaderTest")) {
		// Tu wczytujesz dane z pliku. Możesz też podać ścieżkę jako parametr,
		// wtedy sie nie pokaże dialog wyboru pliku
		// Ważne:
		// 1. zakładam, że dane sa rozdzielone średnikami
		// 2. pierwsza linia musi zawierać ngłówki - czyli nazwy rzeczywistych cech
		// Jeśli miałbyś dane nie spełniające, zwłaszcza tego drugiego, to będziemy kombinować
		CsvReader::CsvData csv_data = CsvReader::loadFile();

		// Kazda linia reprezentuje grupę, pierwszy string to jej nazwa
		// następnie w {} są etykiety zdefiniwanych przez ciebie wymiarów w tej grupie
		// Jesli {} jest puste tzn, że grupa nie moze mieć etykiet (np.: nienazwane)
		// UWAGA - to nie sa nazwy wymiarów związanych z plikiem źródłowym (np. wino)
		// tylko to jak je będziesz używał podczas wizualizowalizacji
		// powiązanie tego z np. winem dzieje sie dopiero w okienku dialogowym
		// dlatego jest to niezależne od rodzaju danych
		// Możesz sobie to dowolnie zmodyfikować, np dodac nowe elementy
		QVector<GroupDefinition> groupDefs = {
			{ "unnamed", {} },
			{ "spacial", { "X", "Y", "Z", "T" } },
			{ "visual", { "Skin_C", "Hair_C", "Eye_S", "Nose_L", "Mouth_W", "Smile", "Frown", "Hair_L", "Face_Elong", "Iris_C"}},
		};
		/* to dopisz: */
		QVector<QPair<int, int>> defs = {
	{0,-1}, // pierwsza cyfra to indeks grupy z listy powyżej (tu np. unnamed=0, spacial=1 itd)
	{1,0},  // a druga cyfra to indeks etykiety (też j.w.), dla grupy unnamed jest pomijany, wpisałem -1
	{1,1},  // tych linii może być mniej niż parametrów win,
	{0,-1},  // wtedy ostatnie pola w dialogu nie bedą po prostu ustawione
	{2,0},  // może być też więcej, wtedy one będą pominięte
	{2,1},
	{2,2},
	{2,3},
	{2,4},
	{0,-1},
	{0,-1}
		};

		defs = {
			/* fix.ac.ty */	{2,9},
			/* vol.ac.ty */	{2,1},
			/* citric a. */	{0,-1},
			/* res.sugar */	{1,2},
			/* chlorides */	{0,-1},
			/* free SO2  */	{2,6},
			/* total SO2 */	{0,-2},
			/* density   */	{2,0},
			/* pH        */ {1,0},
			/* sulphates */	{0,-1},
			/* alcohol   */	{2,3},
			/* quality   */	{2,5},
		};

#if 0
		defs = {
			/* fix.ac.ty */	{,},
			/* vol.ac.ty */	{,},
			/* citric a. */	{,},
			/* res.sugar */	{,},
			/* chlorides */	{,},
			/* free SO2  */	{,},
			/* total SO2 */	{,},
			/* density   */	{,},
			/* pH        */ {,},
			/* sulphates */	{,},
			/* alcohol   */	{,},
			/* quality   */	{,},
		};
#endif

		// Tu tworzysz sobie okienko dialogowe o takich parametrach jak ustawiłeś wyżej
		// csv_data.headers - to sa rzeczywiste etykiety z pliku (np. kwasowośc)
		// groupDefs - to są opcje które będziesz mógł wybierać, ustawione wyżej
		// // STARE:  CsvColumnAssignmentDialog dlg(csv_data.headers, groupDefs);
		CsvColumnAssignmentDialog dlg(csv_data.headers, groupDefs, defs);

		// Tu uruchamiasz dialog i sprawdzasz czy kliknięto OK
		if (dlg.exec() == QDialog::Accepted) {

			// Przypisań i macierzy mozesz teraz używać w dowolnym miejscu plugina
			// nie musisz tego robic jednym ciągiem w tym miejscu
			// radze tylko zawsze na początku własnej procedury sprawdzić czy nie są puste
			// (tzn -> czy plik został wczytany i pogrupowany)
			current_assignment = dlg.getAssignments();




			// Konwersja do macierzy Eigen::Matrixxd, na której możesz sobie dalej działać.
			// Ta macierz ma tyle kolumn co wymiarów, a wierszy tyle co próbek w źródle z pominięciem nagłówków
			current_data_matrix = CsvReader::convertCsvDataToMatrix(csv_data);
			// Ta macierz jest widoczna w całym pluginie więc możesz sobie inne kawałki kodu
			// stworzyć jako reakcję na inny przycisk i ją tam odczytać. Warto tylko sprawdzić czy nie jest pusta.
			ConcretePlugin::K3FormProjectionMatrix(&current_data_matrix);
			// Tu LL: Wypełniam częściowo current_assignment, żeby mieć mniej klikania:
			// current_assignment[0].groupIndex = 1;
			// current_assignment[0].label_id = 0;  // Nie wyszło, odpuszczam.

		}
		else {
			// Jeśli nie kliknieto OK - to zakładam, ze zrezygnowałeś i odrzucam wszystko co wczytałeś
			current_data_matrix = Eigen::MatrixXd();
			current_assignment.clear();
		}

	}
	else if (0 == name.compare(L"K3CzteroPajak")) {

#define h05 0.3536
		// sqrt(2)/4

		CVector3d LiPa(0.0, 0.0, 2.0);
		double scale = 70.0;
		CModel3D* K3MyModel = new CModel3D();
		double QuadStem[3] = { scale * -0.7, scale * 2., scale * 1.8 },
			QuadX[3] = { scale * (-0.5), scale * 0., scale * h05 },
			QuadY[3] = { scale * 0.5, scale * 0, scale * h05 },
			QuadZ[3] = { 0, scale * 0.5, scale * (-h05) },
			QuadT[3] = { 0.0, scale * (-0.5), scale * (-h05) },
			QX[3], QY[3], QZ[3], QT[3];

		CRGBA K3Red((unsigned char)255, (unsigned char)0, (unsigned char)0);
		CRGBA K3Green((unsigned char)0, (unsigned char)255, (unsigned char)0);
		CRGBA K3Blue((unsigned char)100, (unsigned char)100, (unsigned char)255);
		CRGBA K3Yellow((unsigned char)255, (unsigned char)255, (unsigned char)0);

		for (int i = 0; i < 3; i++) {
			QX[i] = QuadX[i];
			QuadX[i] += QuadStem[i];
			QY[i] = QuadY[i];
			QuadY[i] += QuadStem[i];
			QZ[i] = QuadZ[i];
			QuadZ[i] += QuadStem[i];
			QT[i] = QuadT[i];
			QuadT[i] += QuadStem[i];
		};
		double K3R = 10.0;
		K3Arrow* Arro1 = new	K3Arrow(QuadStem, QuadX, K3R, &K3Red);
		K3MyModel->addChild(Arro1);
		K3Arrow* Arro2 = new	K3Arrow(QuadStem, QuadY, K3R, &K3Green);
		K3MyModel->addChild(Arro2);
		K3Arrow* Arro3 = new	K3Arrow(QuadStem, QuadZ, K3R, &K3Blue);
		K3MyModel->addChild(Arro3);
		K3Arrow* Arro4 = new	K3Arrow(QuadStem, QuadT, K3R, &K3Yellow);
		K3MyModel->addChild(Arro4);

		//			0, 1, 2, 3,  // 1 2 swapped
		//			0, 2, -1, 3,  // 1 3
		//			0, 3, -1, -2,  // 0 2DR
		//			-1, -2, -0, -3,  // 0 3
		//			-3, -2, -0, 1
		double NX[3], NY[3], NZ[3], NT[3];
		K3Neg3(QX, NX);
		K3Neg3(QY, NY);
		K3Neg3(QZ, NZ);
		K3Neg3(QT, NT);
		K3R = 4.0;
		for (int jarc = 0; jarc < 5; jarc++) { // emulate rolling
			K3ArrowsArc(QuadStem, NT, NZ, K3MyModel, K3R + .0, 8, &K3Blue);
			K3ArrowsArc(QuadStem, NZ, NY, K3MyModel, K3R + .0, 8, &K3Green);
			K3ArrowsArc(QuadStem, NY, NX, K3MyModel, K3R + .0, 8, &K3Red);
			K3ArrowsArc(QuadStem, NX, NT, K3MyModel, K3R + .0, 8, &K3Yellow);

		};

		K3MyModel->importChildrenGeometry();
		AP::WORKSPACE::addModel(K3MyModel);
	} // End K3CzteroPajak
	else if (0 == name.compare(L"K3Krata"))
	{
		CModel3D* K3MyModel = new CModel3D();
#define Nkrat 8
#define Mkrat 6
		for (int j = 0; j < Mkrat; j++) {
			for (int i = 0; i < Nkrat; i++) {
				K3Totem* ToTen1;
				Eigen::VectorXd K3HyperSpot(20);
				Eigen::VectorXd K3HyperLook(20);
				for (int i88 = 0; i88 < 20; i88++) {
					K3HyperSpot(i88) = 0.0;
					K3HyperLook(i88) = 0.5;
				};
				K3HyperSpot(0) = ((double)i) * 2;
				K3HyperSpot(1) = ((double)j) * 1.5;
				K3HyperSpot(2) = (double)2;
				int FeatureSel[Nkrat] = { 0,1, 2, 3,5,8,10, 19 };
				//  Feature mapping:
	// +0 skin colour (rainbow scale)
	// +1 hair colour (incl. facial hair, if any)
	// +2 - eye size
	// +3 - nose height
	// +4 - mouth width
	// +5 - smile
	// 6 - eybrow frown
	// +7 - hair length
	// +8 - face elongation
	// +9 - iris color
				K3HyperLook(FeatureSel[i]) = (double)j / (double)Mkrat;
				ToTen1 = new K3Totem(K3HyperSpot, K3HyperLook);
				//K3Totem::K3Totem(Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook) {
				// Create a totem whose position and appearance represents data.
				// Assume the given K3HyperSpot = DataPoint*Observer ,
				// i.e. has been projected into the observation space.
				K3MyModel->addChild(ToTen1);
			};

			// }

		}
		AP::WORKSPACE::addModel(K3MyModel);
		UI::updateAllViews();
		// _Thrd_yield();

	}
	else if (0 == name.compare(L"K3Display"))
	{
		CModel3D* K3MyModel = new CModel3D();
		// Eigen::MatrixXd K3FullCloud = current_data_matrix;
		Eigen::MatrixXd K3DenseCloud = K3BoiledData; // *K3_IObs* K3FullCloud;
		//		AP::WORKSPACE::addModel(K3MyModel);
		AP::WORKSPACE::setCurrentModel(-1);
		MatrixXd K3ViewMat(5, 5); // QQ Assuming we only deal with 4 spatial dimensions
		K3FillMat(K3ViewMat, 0.0);
		for (int i = 0; i < 4; i++) { // QQ Same as above: 4D only
			K3ViewMat(i, i) = 50.0;  // QQ Scale! bigbig
		}
		K3ViewMat(4, 4) = 1.0;
		K3AddMyCloud(K3MyModel, K3DenseCloud, K3ViewMat, 2000.3);
		K3MyModel->importChildrenGeometry();

		AP::WORKSPACE::addModel(K3MyModel); // QQ blok Wlknoc
		UI::updateAllViews();
		// _Thrd_yield();


	}
	else if (0 == name.compare(L"K3Dance"))
	{
		CModel3D* K3MyModel = new CModel3D();
		MatrixXd K3ViewMat(5, 5); // QQ Assuming we only deal with 4 spatial dimensions
		// Eigen::MatrixXd K3FullCloud = current_data_matrix;
		Eigen::MatrixXd K3DenseCloud = K3BoiledData; // *K3_IObs* K3FullCloud;
		//		AP::WORKSPACE::addModel(K3MyModel);
		AP::WORKSPACE::setCurrentModel(-1);
		void* k3viewer = UI::CAMERA::currentViewer();

		delete_old_screenshots("Fot_*.png");


		// QQ9 ROLLING THE VOLUME LOOP
		for (int i_plane = 0; i_plane < 5; i_plane++) {  // was 5 and is 12
			for (double alfa = 0.0; alfa < 0.01 + 3.1415926 / 2.0; alfa += (3.1415926 / 24.0)) {

				{
					FILE* plik = _wfopen(DATA_PATH(L"Liczniki.txt").c_str(), L"a");
					if (plik == NULL) {
						UI::MESSAGEBOX::error(L"To bardzo skomPLIKowane");
					}

					fprintf(plik, "LL=%06d_x_%06d_OTO\n", i_plane,
						(int)(100.0 * alfa + 3000.0));
					fclose(plik);
				}

				// Sleep(200);
				// _Thrd_yield();
				// AP::processEvents();
//					AP::WORKSPACE::removeAllModels();
				AP::processEvents();
				// AP::WORKSPACE::removeImage(;
				// K3MyModel = new CModel3D();
				// double alfa = (double)i * 0.25;  // radians!
// void K3_4x4viewN(MatrixXd *V, int k, double alfa)

				K3_4x4viewN(&K3ViewMat, i_plane, alfa);
				K3ViewMat *= K3GRAND_SCALE;   //QQ9 SCALE
				K3ViewMat(4, 4) = 1.0;
				//{
				//	char DescrIter[300];   // QQ w or no w?  wchar_t

				//	sprintf(DescrIter, "Rot_%03d_x_%04d ", i_plane,
				//		(int)(100.0 * alfa + 3000.0));
					// (L"C:/K3/Wielowymiar/MojeWidoki.txt", K3ViewMat, DescrIter);
				//};

				// if (i % 2)
				{
					double RoseStem[3] = { -30, 40, 25 }; // QQ sET STEM LENGTH; was triple 1.
					AP::WORKSPACE::removeAllModels();

					K3MyModel = new CModel3D();
					AP::WORKSPACE::addModel(K3MyModel);
					AP::WORKSPACE::setCurrentModel(-1);

					// K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3DenseCloud, "ViewDense"); // FotFilNam);
					K3AddMyCloud(K3MyModel, K3DenseCloud, K3ViewMat, 5000.4); // 5000000.3);





					// Znakiem QQ7 oznaczam blokady z 23 IX 2024
					CRGBA colorlist[20];
					// CRGBA K3_color(double en, double a)
					// for (int i = 0; i < 15; i++) {
						// colorlist[i] = K3_color((double(i) / 12.0), ((i > 1) ? 0.2 : 0.9));
						// colorlist[i] = K3_color(K3AverIndex(Observer, i) / 15.0, 0.8);
					//};
					colorlist[0] = CRGBA((unsigned char)255, (unsigned char)0, (unsigned char)0);
					colorlist[1] = CRGBA((unsigned char)0, (unsigned char)255, (unsigned char)0);
					colorlist[2] = CRGBA((unsigned char)0, (unsigned char)0, (unsigned char)255);
					colorlist[3] = CRGBA((unsigned char)255, (unsigned char)255, (unsigned char)0);
					K3RoseOfWinds* MyRose = new K3RoseOfWinds(K3ViewMat, RoseStem, colorlist, 20.0);
					// K3MyModel->addChild(MyRose);

					// QQ Fpentli!
					//					K3RoseO1fWinds* IlNome = new K3RoseOfWinds(Observer, RoseStem, colorlist);
					//					K3MyModel->addChild((CMesh*)IlNome); // QQ future orphan

					K3MyModel->importChildrenGeometry();

					UI::updateAllViews(false);

					qInfo() << "checkpoint" << Qt::endl;

					//AP::processEvents();
					//Sleep(40);

					QString K3QST = QString().sprintf("Fot_%03d_x_%04d_PK.png", i_plane, static_cast<int>(100.0 * alfa + 3000.0));

					qInfo() << "Zrzut ekranu do pliku: " << K3QST << " po nazwie";

					QByteArray ba = K3QST.toLocal8Bit();
					const char* c_str2 = ba.data();

					K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3ViewMat, c_str2); // "ViewMat"); // FotFilNam);

				UI:Beep(440.0, 500.0);

					UI::CAMERA::screenshot(PNG_PATH(K3QST), k3viewer); // To teraz będzie tutaj

				};
			};
		};
		//// :end of K3
		// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3ViewMat, "Po Pentlach!!!");

		AP::WORKSPACE::removeAllModels();
		qInfo() << "To już jest koniec..." << Qt::endl;
	};
};


// }