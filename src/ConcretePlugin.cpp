
#include "ConcretePlugin.h"

#include "AP.h"
#include "AppSettings.h"

#include <Eigen/Geometry>
#include <Eigen/Eigenvalues> 

#include "K3Helpers.h"

#include "K3Totem.h"
#include "K3Arrow.h"
#include "K3ChernoffFace.h"
#include "K3RoseOfWinds.h"


using namespace Eigen;

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

//#include "Image.h"

// Now to the job:

#define K3Dimensionality 15

// :koniec dodatku LL


ConcretePlugin::ConcretePlugin(void)
{
	current_data.headers = QStringList();
	current_data.matrix = Eigen::MatrixXd();
	current_assignment = QVector<ColumnAssignment>();
}

ConcretePlugin::~ConcretePlugin(void)
{
}


void ConcretePlugin::onLoad()
{
	UI::PLUGINPANEL::create(m_ID, L"N-Dim View");

	UI::PLUGINPANEL::addButton(m_ID, L"csvReaderTest", L"Load data from .csv/.dat", 0, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"assignGroups", L"Assign data to groups", 1, 0);

	UI::PLUGINPANEL::addButton(m_ID, L"K3Display", L"Static view", 6, 0);
	//UI::PLUGINPANEL::addButton(m_ID, L"K3CzteroPajak", L"Make four-spider", 7, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3Krata", L"Create a grid of avatars", 8, 0);
	UI::PLUGINPANEL::addButton(m_ID, L"K3Dance", L"Let's dance", 9, 0);
}



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


//========================================================
// RawData is probably current_data.matrix
// Result is K3BoiledData
//========================================================
int ConcretePlugin::K3FormProjectionMatrix(const Eigen::MatrixXd &RawData) {
	// Use list of parameter choices and data matrix to create a projection matrix
	// that projects the raw dimensions into the observables according to user wishes
	int K3I_ListAnon[100];
	int i88, j88, j_raw, j;
	int K3I_taken_s[4] = { 0,0,0,0 };
	int K3I_taken_v[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int K3I_Count_a = 0, K3I_Count_s = 0, K3I_Count_v = 0;
	if (RawData.rows() * RawData.cols() > 0)
	{
		int K3_IDim = RawData.rows();
		int K3_ICard = RawData.cols();

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
				X_anon(j88, i88) = RawData(j_raw, i88);
			};
		};
		if (K3I_Count_s < 4) {
			int j_boiled = 0;
			MatrixXd X88 = K3_Get_PCA_Funnel(X_anon, 4 - K3I_Count_s); // QQ K3I_Count_a);
			K3ListMatrix(DATA_PATH("MojeLeje.txt"), X88, "Lejek");
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

		K3BoiledData = *K3_IObs * RawData;
		K3ListMatrix(DATA_PATH("MujStat.txt"), RawData, "RawData");
		K3ListMatrix(DATA_PATH("MujStat.txt"), *K3_IObs, "K3_IObs");
		K3ListMatrix(DATA_PATH("MujStat.txt"), K3BoiledData, "K3BoiledData");
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
	K3ListMatrix(DATA_PATH("MujStat.txt"), *K3_IObs, "NEW_K3_IObs");

	return(1);
};



void saveFileData(QSettings* settings, const QString& filePath, const QVector<QPair<int, int>>& vec)
{
	QString key = "assignments/"+QString::fromUtf8(filePath.toUtf8().toBase64(QByteArray::Base64UrlEncoding));

	settings->remove(key);

	settings->beginGroup(key);

	QVariantList list;
	for (const auto& pair : vec)
		list << QVariant::fromValue(QPoint(pair.first, pair.second)); // łatwo zapisywalne

	settings->setValue("definitions", list);
	settings->setValue("timestamp", QDateTime::currentDateTimeUtc().toString(Qt::ISODate));
	settings->endGroup();

	settings->sync();
}

QVector<QPair<int, int>> loadFileData(QSettings* settings, const QString& filePath)
{
	QVector<QPair<int, int>> result;
	
	QString key = "assignments/" + QString::fromUtf8(filePath.toUtf8().toBase64(QByteArray::Base64UrlEncoding));


	settings->beginGroup(key);
	QVariantList list = settings->value("definitions").toList();
	for (const auto& item : list) {
		QPoint p = item.toPoint();
		result.append({ p.x(), p.y() });
	}
	settings->endGroup();

	return result;
}


void cleanupOldEntries(QSettings* settings, int daysOld = 30)
{
	QDateTime now = QDateTime::currentDateTimeUtc();

	settings->beginGroup("assignments");
	QStringList keys = settings->childGroups();
	for (const QString& group : keys) {
		settings->beginGroup(group);
		QString tsStr = settings->value("timestamp").toString();
		QDateTime ts = QDateTime::fromString(tsStr, Qt::ISODate);
		if (ts.isValid() && ts.daysTo(now) > daysOld)
			settings->remove(""); // usuwa bieżącą grupę
		settings->endGroup();
	}
	settings->endGroup();
	settings->sync();
}


void ConcretePlugin::AssignGroups() {
	if (current_data.matrix.size() == 0) {
		UI::MESSAGEBOX::error("You need to load data first !", "N-Dim-view plugin error");
		return;
	}

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
	
	
	// pierwsza cyfra to indeks grupy z listy powyżej (tu np. unnamed=0, spacial=1 itd)
	// a druga cyfra to indeks etykiety (też j.w.), dla grupy unnamed jest pomijany, wpisałem -1
	// tych linii może być mniej niż parametrów win,
	// wtedy ostatnie pola w dialogu nie bedą po prostu ustawione
	// może być też więcej, wtedy one będą pominięte
	//QVector<QPair<int, int>> defs = {
	//	/* fix.ac.ty */	{2,9},
	//	/* vol.ac.ty */	{2,1},
	//	/* citric a. */	{0,-1},
	//	/* res.sugar */	{1,2},
	//	/* chlorides */	{0,-1},
	//	/* free SO2  */	{2,6},
	//	/* total SO2 */	{0,-2},
	//	/* density   */	{2,0},
	//	/* pH        */ {1,0},
	//	/* sulphates */	{0,-1},
	//	/* alcohol   */	{2,3},
	//	/* quality   */	{2,5},
	//};

	auto settings = AppSettings::pluginSettings("N-Dim-view");
	

	//saveFileData(settings.get(), current_data.file_path, defs);
	QVector<QPair<int, int>> defs = loadFileData(settings.get(), current_data.file_path);

	// Tu tworzysz sobie okienko dialogowe o takich parametrach jak ustawiłeś wyżej
	// csv_data.headers - to sa rzeczywiste etykiety z pliku (np. kwasowośc)
	// groupDefs - to są opcje które będziesz mógł wybierać, ustawione wyżej
	// // STARE:  CsvColumnAssignmentDialog dlg(csv_data.headers, groupDefs);
	CsvColumnAssignmentDialog dlg(current_data.headers, groupDefs, defs);

	// Tu uruchamiasz dialog i sprawdzasz czy kliknięto OK
	if (dlg.exec() == QDialog::Accepted) {

		// Przypisań i macierzy mozesz teraz używać w dowolnym miejscu plugina
		// nie musisz tego robic jednym ciągiem w tym miejscu
		// radze tylko zawsze na początku własnej procedury sprawdzić czy nie są puste
		// (tzn -> czy plik został wczytany i pogrupowany)
		current_assignment = dlg.getAssignments();

		defs.clear();
		for (auto ass : current_assignment) {
			defs.push_back(QPair<int, int>(ass.groupIndex, ass.label_id.value_or(-1)));
		}

		if (!defs.empty()) {
			saveFileData(settings.get(), current_data.file_path, defs);
		}

		K3FormProjectionMatrix(current_data.matrix);
	}
	else {
		// Jeśli nie kliknieto OK - to zakładam, ze zrezygnowałeś i odrzucam wszystko co wczytałeś
		current_data.headers.clear();
		current_data.matrix = Eigen::MatrixXd();
		current_assignment.clear();
	}
}



void ConcretePlugin::K3CzteroPajak(double h05) {

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
}


void ConcretePlugin::K3Krata(int Nkrat, int Mkrat) {
	CModel3D* K3MyModel = new CModel3D();
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

			std::vector<int> FeatureSel( { 0, 1, 2, 3, 5, 8, 10, 19 } ); // Note: size should be equal to Nkrat

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


void ConcretePlugin::K3Display() {
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


void ConcretePlugin::K3Dance(double grand_scale) {
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
				FILE* plik = fopen(DATA_PATH("Liczniki.txt").toStdString().c_str(), "a");
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
			K3ViewMat *= grand_scale;   //QQ9 SCALE
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

				K3ListMatrix(DATA_PATH("MujZrzut.txt"), K3ViewMat, c_str2); // "ViewMat"); // FotFilNam);

				//UI::Beep(440.0, 500.0);

				UI::CAMERA::screenshot(PNG_PATH(K3QST), k3viewer); // To teraz będzie tutaj

			};
		};
	};
	//// :end of K3
	// K3ListMatrix(L"C:/K3/Wielowymiar/MujZrzut.txt", K3ViewMat, "Po Pentlach!!!");

	AP::WORKSPACE::removeAllModels();
	qInfo() << "To już jest koniec..." << Qt::endl;
}


void ConcretePlugin::onButton(std::wstring name)
{
	if (0 == name.compare(L"csvReaderTest")) {
		current_data = CsvReader::loadFile();
		current_assignment = QVector<ColumnAssignment>();
	}
	else if (0 == name.compare(L"assignGroups")) {
		AssignGroups();
	}
	else if (0 == name.compare(L"K3CzteroPajak"))
	{
		K3CzteroPajak(0.3536); //parametr zamiast: #define h05 0.3536
	}
	else if (0 == name.compare(L"K3Krata"))
	{
		// parametry zamiast: #define Nkrat 8 i #define Mkrat 6
		K3Krata(8, 6);
	}
	else if (0 == name.compare(L"K3Display"))
	{
		K3Display();
	}
	else if (0 == name.compare(L"K3Dance"))
	{
		K3Dance(K3GRAND_SCALE);
	}
}

