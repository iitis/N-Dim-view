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

#define K3GRAND_SCALE 5.0
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



//void K3ReadCSV_WithHeader(MatrixXd* X, const wchar_t* filnam) {
//	char K3ParNames[20][30];
//
//	int i;
//	int j = 0;
//	int nmin = 20;
//
//	for (i = 0; i < 20; i++) {
//		K3ParNames[i][0] = (char)0;
//	}
//	FILE* plik = _wfopen(filnam, L"r");
//	if (plik) {
//		int DoneRead = 8;
//		j = 0;
//		MatrixXd ThisSpot(1, 20);
//		/* DoneRead=fscanf(plik, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s",
//		K3ParNames[0], K3ParNames[1], K3ParNames[2], K3ParNames[3], K3ParNames[4],
//		K3ParNames[5], K3ParNames[6], K3ParNames[7], K3ParNames[8], K3ParNames[9],
//		K3ParNames[10], K3ParNames[11], K3ParNames[12], K3ParNames[13], K3ParNames[14],
//		K3ParNames[15], K3ParNames[16], K3ParNames[17], K3ParNames[18], K3ParNames[19]); */
//		while (DoneRead != EOF) {
//			char K3buf[300];
//			char* K3Yad;
//			K3Yad = fgets(K3buf, 290, plik);
//			if (K3Yad == K3buf) {
//				DoneRead = sscanf(K3buf, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf",
//					&ThisSpot(0, 0), &ThisSpot(0, 1), &ThisSpot(0, 2), &ThisSpot(0, 3), &ThisSpot(0, 4),
//					&ThisSpot(0, 5), &ThisSpot(0, 6), &ThisSpot(0, 7), &ThisSpot(0, 8), &ThisSpot(0, 9),
//					&ThisSpot(0, 10), &ThisSpot(0, 11), &ThisSpot(0, 12), &ThisSpot(0, 13), &ThisSpot(0, 14),
//					&ThisSpot(0, 15), &ThisSpot(0, 16), &ThisSpot(0, 17), &ThisSpot(0, 18), &ThisSpot(0, 19));
//				if (DoneRead != EOF) {
//					K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", ThisSpot, "ThisSpot read");
//					if ((DoneRead > 2) /* && (j < 1400) */) {  // QQ ogran.
//						for (i = 0; i < DoneRead; i++) {
//							(*X)(i, j) = ThisSpot(0, i);
//						}
//						if (DoneRead < nmin) {
//							nmin = DoneRead;
//						}
//					}
//					else
//					{
//						DoneRead = EOF;
//					}
//					j++;
//				}
//			}
//			else
//			{
//				DoneRead = EOF;
//			}
//		}
//	};
//	MatrixXd A;
//	A = X->block(0, 0, nmin, j);
//	*X = A; // .transpose();
//
//	fclose(plik);
//};


//void K3LogEntry(const wchar_t* filnam, const char* mujText) {
//	int i, j;
//	time_t czas;
//	FILE* plik = _wfopen(filnam, L"a");
//	if (plik == NULL) {
//		UI::MESSAGEBOX::error(L"To bardzo skomPLIKowane");
//	}
//	time(&czas);
//	fprintf(plik, "LOG at Date=%s : %s \n", ctime(&czas), mujText);
//	fclose(plik);
//};

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

	UI::PLUGINPANEL::addButton(m_ID, L"K3MetaChmura", L"Zrób K3MetaChmurę", 6, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3CzteroPajak", L"Zrób K3CzteroPajak", 7, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3Krata", L"Zrób Kratę", 8, 0);
}


//#include "AnnotationPlane.h"
//#include "AnnotationPoint.h"
//#include <vector>
//#include <Eigen/src/Eigenvalues/EigenSolver.h>



//void K3FillMeshToUnitCube(CMesh* ThisMesh, CRGBA Kolor) { // double l, double a) {
//	CMesh* K3UnitCube = new CMesh();
//	ThisMesh->addVertex(CVertex(0.0, 0.0, 0.0)); // , Kolor);
//	ThisMesh->addVertex(CVertex(0.0, 0.0, 1.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(0.0, 1.0, 0.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(0.0, 1.0, 1.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(1.0, 0.0, 0.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(1.0, 0.0, 1.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(1.0, 1.0, 0.0)); //, Kolor);
//	ThisMesh->addVertex(CVertex(1.0, 1.0, 1.0)); //, Kolor);
//	// looking from outside, vertices should be counterclockwise:
//	ThisMesh->faces().push_back(CFace(1, 5, 3));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(7, 3, 5)); // face 0
//
//	ThisMesh->faces().push_back(CFace(4, 6, 5));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(7, 5, 6)); // face 1
//
//	ThisMesh->faces().push_back(CFace(0, 2, 4));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(6, 4, 2)); // face 2
//
//	ThisMesh->faces().push_back(CFace(1, 3, 0));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(2, 0, 3)); // face 3
//
//	ThisMesh->faces().push_back(CFace(2, 3, 6));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(7, 6, 3)); // face 4
//
//	ThisMesh->faces().push_back(CFace(0, 4, 1));
//	//ThisMesh->fcolors().push_back(Kolor);
//	ThisMesh->faces().push_back(CFace(5, 1, 4)); // face 5
//
//	// UI::MESSAGEBOX::error(L"Gonna color");
//
//	// Make front face translucent and red:
//	for (int i = 0; i < 2; i++) {
//		// ThisMesh->fcolors().push_back(CRGBA(0.9f, .3f, .3f, 0.3f));
//		ThisMesh->fcolors().push_back(Kolor);
//		// ThisMesh->fcolors().push_back(K3Color);
//	};
//
//	for (int i = 0; i < 10; i++) {
//		// ThisMesh->fcolors().push_back(CRGBA(0.3f, 0.3f, 1.0f, 1.0f));
//		ThisMesh->fcolors().push_back(Kolor);
//		// ThisMesh->fcolors().push_back(K3Color);
//	};
//
//};

// {
//	CImage* im = new CImage(cx, cy, CImage::Format::Format_ARGB32);

//
//void K3AddTotemMeshFromSpotQQQnicWeWogle(CMesh* whereto, Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook) {  // "HyperSight = DataPoint*Observer"
//	CMesh* K3Korpus = new CMesh;
//	Eigen::Matrix4d K3Position;
//	K3Position << 1, 0, 0, K3HyperSpot[0],
//		0, 2, 0, K3HyperSpot[1],
//		0, 0, .6, K3HyperSpot[2],
//		0, 0, 0, 1;
//
//	/*Eigen::Matrix4d K3MakeKorpus;
//	K3MakeKorpus << 1, 0, 0, 0,
//		0, 2, 0, 0,
//		0, 0, .2, 0,
//		0, 0, 0, 1;*/
//
//	K3FillMeshToUnitCube(K3Korpus, CRGBA(1.0f, 1.0f, 0.70f, 0.990f));  // qq kolorki
//
//	// Scale up the torso:
//	for (int i = 0; i < K3Korpus->vertices().size(); i++)
//	{
//		K3Korpus->vertices()[i].transformByMatrix(K3Position);
//	}
//
//
//	for (int i = 0; i < K3Korpus->vertices().size(); i++)
//	{
//		K3Korpus->vertices()[i].transformByMatrix(K3Position);
//	}
//
//	for (int i = 0; i < K3Korpus->vertices().size(); i++) {
//		whereto->vertices().push_back(K3Korpus->vertices()[i]);
//	};
//
//	for (int i = 0; i < K3Korpus->faces().size(); i++) {
//		whereto->faces().push_back(K3Korpus->faces()[i]);
//	};
//
//	// delete korpus;
//
//
//	Eigen::Matrix4d K3MakeArm;
//	K3MakeArm << 5, 0, 0, 0,
//		0, .2, 0, 0,
//		0, 0, .05, 0,
//		0, 0, 0, 1;
//
//	double alfa = K3HyperLook[4];
//
//	Eigen::Matrix4d K3TwistArm;
//	K3TwistArm << cos(alfa), -sin(alfa), 0, 0,
//		sin(alfa), cos(alfa), 0, 0,
//		0, 0, 1, 0,
//		0, 0, 0, 1;
//
//	Eigen::Matrix4d K3ShiftArm;
//	K3ShiftArm << 1, 0, 0, 2,
//		0, 1, 0, 2,
//		0, 0, 1, 0,
//		0, 0, 0, 1;
//
//};
//

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

void  ConcretePlugin::K3MetaChmura(CModel3D* o, float d)
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
	int m = K3ObsCloud.rows();
	MatrixXd K3LocalSpots(K3hs + K3asr + 1, k);  // WRONG: 5 = X Y Z d T (3D + distance from current 3D space + homogeneous extra dim)
			// RIGHT??: 5 = X Y Z T d (4D + distance from current 3D space)
	for (int i = 0; i < k; i += 1) {  // DARKU TUTAJ DALEM 10, ZEBY WYSWIETLAC MNIEJ ELEMENTOW
		for (int j = 0; j < K3hs + K3asr + 1; j++) {
			// std::cout << "shame on i, j, k:" << i << " " << j << " " << k << std::endl;
			if (j < m) {  // missing dimensions? Fill with ones.
				K3LocalSpots(j, i) = K3ObsCloud(j, i);
			}
			else {
				K3LocalSpots(j, i) = 1.0;
			}
		}
		K3LocalSpots(K3hs + K3asr, i) = 1.0;  // homogeneous extra dim  // K3UNBLOCKED 2024.08.22
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
		if ((abs(K3LocalSpots(3, i)) < K3Toler) && (abs(K3LocalSpots(4, i)) != 0)) {
			K3Totem* K3TenTotem;
			i = i; //qq
			Eigen::VectorXd P(4);
			Eigen::VectorXd V(K3hv + K3avr);
			// QQSize P is the position vector!
			double askala = 0.1;
			P << K3LocalSpots(0, i) * askala, K3LocalSpots(1, i)* askala, K3LocalSpots(2, i)* askala, K3LocalSpots(3, i)* askala; // Totem
			// P << 15.0 * i, 15.0 * i + 0.1, 15.0 * i + 0.2, 4.04; // Totem
			for (int j88 = 0; j88 < K3hv + K3avr; j88++) {
				V(j88) = K3ObsCloud(5 + j88, i);
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

	K3ListMatrix(DATA_PATH(L"M1120").c_str(), K3ObsCloud, "CloudAdded");
}

void K3Neg3(double A[3], double B[3]) {
	for (int i = 0; i < 3; i++) {
		B[i] = -A[i];
	}
}

void ConcretePlugin::onButton(std::wstring name)
{
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
			{ "spacial", { "dimension_1", "dimension_2", "dimension_3", "dimension_4" } },
			{ "visual", { "left eye", "right eye", "nose", "bgcolor"}},
		};


		QVector<QPair<int, int>> defs = {
			{0,-1}, // pierwsza cyfra to indeks grupy z listy powyżej (tu np. unnamed=0, spacial=1 itd) 
			{1,0},  // a druga cyfra to indeks etykiety (też j.w.), dla grupy unnamed jest pomijany, wpisałem -1
			{1,1},  // tych linii może być mniej niż parametrów win,
			{1,2},  // wtedy ostatnie pola w dialogu nie bedą po prostu ustawione
			{2,0},  // może być też więcej, wtedy one będą pominięte
			{2,1},
			{2,2},
			{2,3},
			{1,3},
			{0,-1},
			{0,-1}
		};

		// Tu tworzysz sobie okienko dialogowe o takich parametrach jak ustawiłeś wyżej
		// csv_data.headers - to sa rzeczywiste etykiety z pliku (np. kwasowośc)
		// groupDefs - to są opcje które będziesz mógł wybierać, ustawione wyżej
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
	else if (0 == name.compare(L"K3MetaChmura"))
	{
		// K3 My Menu Item:
#define K3Dim 15
// #define K3CloudCard 15
#define K3CloudCard 12000
		// 5662 

		delete_old_screenshots("Fot_*.png");

		// Eigen::MatrixXd K3FullCloud(K3Dim, K3CloudCard);
		Eigen::MatrixXd K3FullCloud = current_data_matrix;
		int i88, j88, j_raw;
		int K3I_taken_s[4] = { 0,0,0,0 };
		int K3I_taken_v[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int K3I_Count_a = 0;
		int K3I_ListAnon[100];
		for (i88 = 0; i88 < 100; i88++) {
			K3I_ListAnon[i88] = -1;
		}
		if (K3FullCloud.rows() * K3FullCloud.cols() > 0) {
			int K3_IDim = K3FullCloud.rows();
			int K3_ICard= K3FullCloud.cols();

			Eigen::MatrixXd K3_IObs(4 + 10, K3_IDim); // 4 spatial + 10 face params.
			K3FillMat(K3_IObs, 0.0);

			for (int j = 0; j < K3_IDim; j++) {
				int j_raw= current_assignment[j].featureIndex;
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
					}
					K3_IObs(j_boiled, current_assignment[j].featureIndex) = 1;
					break;
				case 2:  // visual
					if (current_assignment[j].label_id.has_value()) {
						j_boiled = current_assignment[j].label_id.value() + 4;
					};
					K3_IObs(j_boiled, current_assignment[j].featureIndex) = 1;
					break;

				}
			};
			// make funnel. Firt, get the anonymous together:
			MatrixXd X_anon(K3I_Count_a, K3_ICard);
			for (j88 = 0; j88 < K3I_Count_a; j88++) {
				j_raw = K3I_ListAnon[j88];
				for (i88 = 0; i88 < K3_ICard; i88++) {
					X_anon(j88, i88) = K3FullCloud(j_raw, i88);
				};
			};
			MatrixXd X88 = K3_Get_PCA_Funnel(X_anon, K3I_Count_a);



			MatrixXd Observer(K3Dim, K3Dim); // Observer pojdzie precz, zamiast tego cztery macierze.
			// K3Dim shall be split into 4 parts, not necessarily disjoint. See article.

			/*
			// number of human-identifiable parameters to be displayed hyperspacially:
	#define K3hs 2  // density, pH
	int K3_indhs[k3hs] = {7, 8};  // indices in input matrix
	// number of human-identifiable parameters to be displayed by avatar appearance:
	#define K3hv 2  // citric, suger
	int K3_indhv[k3hv] = {2,3};
	// number of NOT human-identifiable parameters to be displayed hyperspacially:
	#define K3as 3
	int K3_indas[k3hv] = {2,3};
	// reduces to:
	#define K3asr 2
	// number of NOT human-identifiable parameters to be displayed by avatar appearance:
	#define K3av 6
	// reduces to:
	#define K3avr 4 */

#define K3hs 2  // density, pH
#define K3as 3
#define K3asr 2
#define K3hv 2  // citric, suger,
#define K3av 6
#define K3avr 4

			MatrixXd K3DenseCloud(K3hs + K3asr + K3hv + K3avr, K3CloudCard); // Wlknoc przeniesiono z dołu
			MatrixXd Obser_s(K3hs + K3asr, K3Dim);
			K3FillMat(Obser_s, 0.0);
			MatrixXd Obser_v(K3hv + K3avr, K3Dim);
			K3FillMat(Obser_v, 0.0);
			Obser_s(0, 7) = 1;
			Obser_s(1, 8) = 1;

			Obser_v(0, 2) = 1;
			Obser_v(1, 3) = 1;


			Obser_s(0, 0) = 1;
			Obser_s(1, 1) = 1;
			Obser_s(2, 4) = 1;
			Obser_s(3, 5) = 1;

			Obser_v(0, 6) = 1;
			Obser_v(1, 9) = 1;
			Obser_v(2, 10) = 1;
			Obser_v(3, 11) = 1;
			//		MatrixXd K3_MessS(K3as, )


			/*		int j = 0;
					for (int i = 0; i < K3asr; i++) {
						int k3flag = 0;
						repeat{ k3flag = 0;
						for(i88;i88<K3hs+K3a)
						}until(2 > 3);
					} */

			int np_hs = 0, np_hv = 0, np_as = 0, np_av = 0;
			// #define K3addHS(kk)  Obser_hs(np_hs++,kk)+=1.0;
					/*	0	fixed_acidity	*/
					/*	1	volatile_acidity	*/
				/*	2	citric_acid	*/
				/*	3	residual_sugar	*/
					/*	4	chlorides	*/
					/*	5	free_sulfur_dioxide	*/
					/*	6	total_sulfur_dioxide	*/
				/*	7	density	*/
				/*	8	pH	*/
					/*	9	sulphates	*/
					/*	10	alcohol	*/
					/*	11	quality	*/




			{
				// K3 loading CSV file:
				char K3ParNames[12][30];


				{// Clear log files: QQQ
					FILE* plik = _wfopen(DATA_PATH(L"K3av.txt").c_str(), L"w");
					fclose(plik);
					plik = _wfopen(DATA_PATH(L"MujZrzut.txt").c_str(), L"w");
					fclose(plik);
					plik = _wfopen(DATA_PATH(L"MojeWidoki.txt").c_str(), L"w");
					fclose(plik);

				}

				// K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3FullCloud, DATA_PATH(L"CircleSquareEtc.dat").c_str());
				// K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3FullCloud, DATA_PATH(L"Tesseract.dat").c_str());
				// K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3FullCloud, "CSVK3FullCloud");
				// UI::MESSAGEBOX::error(L"Got may data");  // LABEL GOTDATA


				// Eigen::VectorXd X(15)[15]; //  double X[5][15];
				std::vector< Eigen::VectorXd> X;
				int j;
				char K3Bufor[80];
				Eigen::Matrix3d m3;
				m3(1, 2) = 742;
				QString K3QST;

				// Test PCA:

				if (1 < 0) {  // parentheses to block range of variables
					MatrixXd X77(3, 100);
					MatrixXd K3Axes77(3, 3);
					K3Axes77 << 3.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 0.0, 5.0;
					for (int i77 = 0; i77 < 100; i77++) {
						double a, b, c;
						a = ((double)(rand() % 10000)) / 10000.0 * 3;
						b = ((double)(rand() % 10000)) / 10000.0 * .8;
						c = ((double)(rand() % 10000)) / 10000.0 * .5;

						X77(0, i77) = a * K3Axes77(0, 0) + b * K3Axes77(0, 1) + c * K3Axes77(0, 2);
						X77(1, i77) = a * K3Axes77(1, 0) + b * K3Axes77(1, 1) + c * K3Axes77(1, 2);
						X77(2, i77) = a * K3Axes77(2, 0) + b * K3Axes77(2, 1) + c * K3Axes77(2, 2);


					};
					FILE* plik = _wfopen(DATA_PATH(L"Gugu.txt").c_str(), L"w");
					fprintf(plik, "X77=[ ");
					for (int i = 0; i < 100; i++) {
						fprintf(plik, "%15lg,  %15lg,  %15lg\n", X77(0, i), X77(1, i), X77(2, i));
					}
					fprintf(plik, "]\n plot3(X77(:,1),X77(:,2),X77(:,3))\n");
					fclose(plik);
					//				std::cout << "X77:" << std::endl << X77 << "X77 over" << std::endl;
					//MatrixXd ConcretePlugin::K3_Get_PCA_Funnel(MatrixXd X, int nd)
					K3ListMatrix(DATA_PATH(L"K3av.txt").c_str(), X77, "L77=1317; X77");
					MatrixXd X88;
					X88 = K3_Get_PCA_Funnel(X77, 2);
					K3ListMatrix(DATA_PATH(L"K3av.txt").c_str(), X88, "L77=1320; X88");
				};  // End of if 1<0





				// FILE* plik = _wfopen(L"C:/K3/Wielowymiar/ala.txt", L"r");   // ok 1599 rekordów
				// FILE* plik = _wfopen(DATA_PATH(L"NewCircleSquareEtc2_17.dat").c_str(), L"r");
				// ok  rekordów

#if 0  // cutting out "inline" data read
				FILE* plik = _wfopen(DATA_PATH(L"ScaledCyrcleSquareEtc2_17.dat").c_str(), L"r");   // ok  rekordów

							// QQQQ Tu naprawde czytamy
				if (plik) {
					// fscanf(plik, "%[^,] %[^,]   %[^,] %[^,] %[^,] %[^,] %[^,] %[^,] %[^,] %[^,] %[^,] %[^,] %lf",
					//fscanf(plik, "%s %s %s %s %s %s %s %s %s %s %s %s",
					// K3ParNames[0], K3ParNames[1], K3ParNames[2],
					//	K3ParNames[3], K3ParNames[4], K3ParNames[5],
					//	K3ParNames[6], K3ParNames[7], K3ParNames[8],
					//	K3ParNames[9], K3ParNames[10], K3ParNames[11]);
					//  fgets(K3Bufor, 80, plik);
					for (int j88 = 0; j88 < K3Dim; j88++) {
						for (int i88 = 0; i88 < K3Dim; i88++) {
							Observer(j88, i88) = 0.0;
						}
						Observer(j88, j88) = 1.0;

					}
					// K3AdjustShare(Observer, int iBigDim, int iSmallDim, double delta);
					K3AdjustShare(Observer, 4, 0, 1000.0);
					K3AdjustShare(Observer, 4, 1, 10.0);
					K3AdjustShare(Observer, 4, 2, 10.0);

					/*
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
	*/
	// Dimension of
#define K3ObsDim (K3hs+K3hv+K3asr+K3avr)

					Eigen::VectorXd K3ThisSpot(K3Dim);
					K3ThisSpot << 30, 2, 1.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

					Eigen::MatrixXd NewObserver(K3ObsDim, K3Dim);

					// for (j = 0; j < 1599; j++) {
					for (j = 0; j < K3CloudCard; j++) {
						// for (int i898 = 0; i898 < 1; i898++) {
						fscanf(plik, "%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf",
							&K3ThisSpot(0), &K3ThisSpot(1), &K3ThisSpot(2), &K3ThisSpot(3), &K3ThisSpot(4),
							&K3ThisSpot(5), &K3ThisSpot(6), &K3ThisSpot(7), &K3ThisSpot(8), &K3ThisSpot(9),
							&K3ThisSpot(10), &K3ThisSpot(11), &K3ThisSpot(12), &K3ThisSpot(13), &K3ThisSpot(14));
						// };
						for (int i = 0; i < K3Dim; i++) {
							K3FullCloud(i, j) = K3ThisSpot(i);
						};
						/*K3ThisSpot(0) *= 10;
						K3ThisSpot(1) *= 10;
						K3ThisSpot(2) *= 10;*/

						VectorXd* X88 = new VectorXd(K3Dim);
						*X88 = Observer * K3ThisSpot;  // QQ 22 II 2024: dodałem * na początku    ######### TRANS #########
						X.push_back(*X88);



						/*X[j][0], X[j][1], X[j][2], X[j][3],
						X[j][4], X[j][5], X[j][6], X[j][7],
						X[j][8], X[j][9], X[j][10], X[j][11]);*/

					}
					fclose(plik);  // UI::MESSAGEBOX::error(L"Got data");
				};
				// :End K3 loading CSV file
#endif
				Eigen::MatrixXd K3GasS(K3as, K3CloudCard);
				// K3ListMatrix(L"C:/K3/Wielowymiar/K3av.txt", K3av, "K3av");
				Eigen::MatrixXd K3GasV(K3av, K3CloudCard);
				K3GasS = Obser_s * K3FullCloud;
				K3ListMatrix(DATA_PATH(L"K3GasS.txt").c_str(), K3GasS, "K3GasS");
				K3GasV = Obser_v * K3FullCloud;
				MatrixXd K3FunS, K3FunV;
				// K3FunS = K3_Get_PCA_Funnel(K3GasS, K3asr);
				// K3FunV = K3_Get_PCA_Funnel(K3GasV, K3avr); //  qq
				K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3GasS, "K3GasS");
				K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3FunS, "K3FunS");
				K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3GasV, "K3GasVS");
				K3ListMatrix(DATA_PATH(L"MujZrzut.txt").c_str(), K3FunV, "K3FunV");

				MatrixXd K3NewObserver(K3hs + K3asr + K3hv + K3avr, K3Dim); // QQ Nie Card!
				K3NewObserver << Obser_s,    // QQ Tutaj byk!
					// K3FunS* Obser_as,
					Obser_v,
					// K3FunV* Obser_av;
				/*K3FillMat(K3NewObserver, 0.0);
				for (int j = 0; j < K3hs; j++) {
					for (int i = 0; i < K3Dim) {
						K3NewObserver(j, i) = Obser_hs(j, i);
					}
				}
				for (int j = K3hs; j < K3hs+ K3asr; j++) {
					for (int i = 0; i < K3Dim) {
						K3NewObserver(j, i) = (j, i);
					}
				}  */

				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", Obser_hs, "Obser_hs");
				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3FunS, "K3FunS");
				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", Obser_hv, "Obser_hv");
				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3FunV, "K3FunV");


				// MatrixXd K3DenseCloud(K3hs + K3asr + K3hv + K3avr, K3CloudCard); // Wlknoc przerzucam na początek
				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3NewObserver, "K3NewObserver");
					K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3FullCloud, "K3FullCloud");
				K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3NewObserver, "K3NewObserver");


				// NEUTRALIZE OBSERVER!
				K3FillMat(K3NewObserver, 0.0);
				for (int j = 0; j < 5; j++) {
					K3NewObserver(j, j) = 1.0;
				}

				K3DenseCloud = K3NewObserver * K3FullCloud;


				// shortcut! QQQQ QQQ QQ 
							// MatrixXd X77(1500, 6);


				// (DP): Zwróć uwagę, że kilka linijek wyżej zapisujesz do tej macierzy wynik mnożenia,
				// wydaje mi się, że któraś z tych operacji jest niepotrzebna.
				// K3DenseCloud = MatrixXd(6, K3CloudCard);
				// (DP): UWAGA NA DRUGI PARAMETR - JEŚLI JEST ZA MAŁY TO SIE WYSYPUJE PRZY WCZYTYWANIU DANYCH


				K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3DenseCloud, L"C:/K3/Wielowymiar/Tabulka.txt");
				// K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3DenseCloud, L"C:/K3/Wielowymiar/NewCircleSquareEtc3_13.dat");


				// (DP): Jak się ma poniższa funkcja do tego że wyżej czytasz z tego samego pliku
				// Czy czasem nie czytasz go dwa razy? NewZigZagEDGE
				// K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3DenseCloud, DATA_PATH(L"NewCyrcleSquareEtc2_17.dat").c_str());
				// K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3DenseCloud, DATA_PATH(L"NewZigZag2_3.dat").c_str());
				// QQQQ K3ReadCSV_WithHeader(/* READING_CLOUD */ &K3DenseCloud, DATA_PATH(L"ScaledCyrcleSquareEtc.dat").c_str());
				// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3DenseCloud, "K3DenseCloud csv");

				MatrixXd K3ViewMat(K3hs + K3asr + 1, K3hs + K3asr + 1);
				K3FillMat(K3ViewMat, 0.0);
				for (int i = 0; i < K3hs + K3asr + 1; i++) {
					K3ViewMat(i, i) = 1.0 * K3GRAND_SCALE;  // QQ Scale!
				}
				// K3ViewMat(K3hs + K3asr, K3hs + K3asr)=1.0;

					//// Trying totem class:  // Wlknoc duzy blok
				CModel3D* K3MyModel = new CModel3D();
				//for (int iDataPoint = 0; iDataPoint < 0; iDataPoint++) { // QQ liczba!
				//	K3Totem* K3ThisTotem = new K3Totem(X[iDataPoint], X[iDataPoint]); //  (K3ThisSpot);
				//	//CMesh K3ThisTotem;
				//	// K3AddTotemMeshFromSpot((CMesh*)K3ThisTotem, K3ThisSpot);

				//	// K3MyModel->addChild((CMesh*)K3ThisTotem); //  .CMesh); // Wlknoc blok
				//	//// delete &K3ThisTotem;
				//	// QQ K3MyModel->importChildrenGeometry();
				//};

				// K3MyModel->importChildrenGeometry();


				//for (int j = 0; j < 3; j++) {  // QQ 777 Debugging trick!
				//	K3DenseCloud(j, 0) = 1.0;
				//	K3DenseCloud(j, 1) = 10.0;
				//	K3DenseCloud(j, 2) = 100.0;
				//}
				//for (int j = 3; j < 4; j++) {  // QQ 777 Debugging trick!
				//	K3DenseCloud(j, 0) = 1.0;
				//	K3DenseCloud(j, 1) = 1.0;
				//	K3DenseCloud(j, 2) = 1.0;
				//}
	//			// (L"C:/K3/Wielowymiar/MujZrzut.txt", K3DenseCloud, "K3DenseCloud");
	//			// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3ViewMat, "K3ViewMat");


				K3AddMyCloud(K3MyModel, K3DenseCloud, K3ViewMat, 2000.3);

				//double RoseStem[3] = { 1, 1, 1 };
				//CRGBA colorlist[20];
				//// CRGBA K3_color(double en, double a)
				//for (int i = 0; i < 15; i++) {
				//	// colorlist[i] = K3_color((double(i) / 12.0), ((i > 1) ? 0.2 : 0.9));
				//	colorlist[i] = K3_color(K3AverIndex(Observer, i) / 15.0, 0.8);
				//};
				// K3RoseOfWinds* IlNome = new K3RoseOfWinds(Observer, RoseStem, colorlist);
				// K3MyModel->addChild((CMesh*)IlNome);
	//			K3MyModel->importChildrenGeometry();
	//			AP::WORKSPACE::addModel(K3MyModel); // QQ blok Wlknoc
	//			UI::updateAllViews();
				// _Thrd_yield();
				// K3MyModel = new CModel3D();

				// GLViewer* k3viewer = UI::CAMERA::currentViewer();
				void* k3viewer = UI::CAMERA::currentViewer();

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
			};

			AP::WORKSPACE::removeAllModels();
			qInfo() << "To już jest koniec..." << Qt::endl;
		}

	}

}