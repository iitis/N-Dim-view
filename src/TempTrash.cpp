/*
MISJA:
Stworzyæ nawigacjê nD

- Za³o¿enia: nie korzystam z ruszania myszk¹, wszystko programowo.

- czy mo¿na reagowaæ na klawisze?

- modyfikacja macierzy po³o¿enia obserwatora:

  -  sterowanie 3D w aktualnej przestrzeni obserwacji;
  - jak siê obracaæ ku innym wymiarom?
	 - do kolumny "ku" dodawaæ porcjami,
		   w wierszach x? y? z?
  - Jak wyœwietlaæ n-wymiarow¹ ró¿ê wiatrów?
	- strza³ki z wersorów tych wymiarów, które maj¹ udzia³ nie mniejszy ni¿...

  -  skrypty w plikach tekstowych;

  CEL NAJBLI¯SZY:
  Trajektoria obserwatora, za³o¿enia na pozycjonowanie w nD.

  - 3 widzialne wymiary ustalaæ w kolejnoœci zgodnej z numeracj¹ w hiperprzestrzeni
  - A je¿eli s¹ mieszane? W kolejnoœci malej¹cej œredniej wa¿onej indeksów w nD.

  Zak³adam - przynajmniej na razie - ¿e wymiary odzwierciedlane wygl¹dem
	  (a nie lokalizaj¹) totemu bêd¹ "pojedyncze", tzn ka¿da cecha wygl¹du
	  oddaje tylko jeden parametr Ÿród³owy. To wa¿ne ograniczenie na macierz Observer.

  Ruchy 3 widzialnych wymiarów w przestrzeni nD:
  PILOTEM nazywam Ÿród³o komend dla obserwatora, mo¿e nim byæ:
	  u¿ytkownik przy klawiaturze,
	  skrypt wczytywany z pliku albo z tablicy sta³ych
  Komend¹ mo¿e byæ:
	  identyfikator (numer? nazwa?) któregoœ wymiaru z nD;
		  tak¹ komendê wykonuje sie przez zwiêkszenie udzia³u wskazanego wymiaru...
			- ale w czym i kosztem czego? Pasowa³oby w tym z Trzech,
			   który ju¿ _tego_ wymiaru ma najwiêcej. Jeœli ¿aden nie ma, to w tym,
			   który jest mu najbli¿ej œredni¹ wa¿on¹ indeksów.
			   I kosztem dotychczasowych sk³adowych tego¿ wymiaru.
			   Tylko: solidarnie czy najchudszy traci?
	  identyfikator (numer? nazwa?) któregoœ wymiaru z 3D;
		  tak¹ komendê wykonuje sie przez obrót widzialnego przekroju
			 wokó³ wskazanego wymiaru o 45 stopni


   PROGRAM PRACY:
   + napisaæ funkcjê, wyznaczaj¹c¹ œredni (wa¿enie) indeks w nD danego wymiaru 3D:
   K3AverIndex(Observer, iSmallDim);
   + Funkcja reguluj¹ca udzia³ danego (indeksem) wymiaru nD w danym (indeksem) w. 3D
   + K3AdjustShare



*/










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
	for (i = 0; i < Observer.cols(); i++) { // QQ mo¿e warto pomin¹æ iBigDim?
		SumOfAbsolutes += abs(Observer(iSmallDim, i));
	};
	if (SumOfAbsolutes > 0.001) {
		Factor = (SumOfAbsolutes - abs(delta)) / SumOfAbsolutes;
		if (Factor < 0.0) {
			Factor = 0.0;
		}
		for (i = 0; i < Observer.cols(); i++) {  // QQ i tutaj te¿?
			Observer(iSmallDim, i) *= Factor;
		};

	};
	Observer(iSmallDim, iBigDim) += delta;
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


