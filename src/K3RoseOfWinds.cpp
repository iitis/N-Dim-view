#include "K3RoseOfWinds.h"

#include "K3Helpers.h"

using namespace Eigen;


// definicja z K3Arrow
#define K3ArrowHeadWidth 0.15


// Now to the job:

#define K3Dimensionality 15
/*K3RoseOfWinds::K3RoseOfWinds(double* RoseStem, CRGBA* colorlist)
{
	double tetra[3][4] = {
		0, 400, 0, 0,
		0, 0, 400, 0,
		0, 0, 0, 400 };

	for (int i = 0; i < 4; i++) {
		this->addVertex(CVertex(tetra[0][i], tetra[1][i], tetra[2][i]));
	};

	this->faces().push_back(CFace(0, 1, 2));
	this->faces().push_back(CFace(0, 3, 1));
	this->faces().push_back(CFace(0, 2, 3));
	this->faces().push_back(CFace(1, 3, 2));

	for (int i = 0; i < 4; i++) {  // looping over the faces:

		this->fcolors().push_back(colorlist[i]);
	};

	MatrixXd Observer777qqq(3,3);
};
*/


K3RoseOfWinds::K3RoseOfWinds(MatrixXd Observer, double* RoseStem, CRGBA* colorlist)
{
	// int j, j88, i;


// CMesh* this = new CMesh();

#define K3ArrowHeadShare 0.1
#define HeadWidth 0.15
#define K3ArrowStemWidth 0.1
	double a = 1 - K3ArrowHeadShare; // the X at neck
	double b = K3ArrowStemWidth;
	double c = K3ArrowHeadWidth;
	double K3UnitArrow[3][10] = {
	0, a, a, a, a, a, a, a, a, 1,
	0, b, b, -b, -b, c, c, -c, -c, 0,
	0, -b, b, b, -b, -c, c, c, -c, 0 };

#define K3VersorScale 200.0  // QQ8 było 200.0. Zmiana na 400 nic dobrego nie zrobila.

	double Xlen; // to check if vector visible


	double BigX[K3Dimensionality], ViewX[K3Dimensionality], SmallX[3], tip[3];


	for (int j = 0; j < K3Dimensionality; j++) {
		for (int j88 = 0; j88 < K3Dimensionality; j88++) {
			BigX[j88] = 0.0;
		};
		BigX[j] = 1.0;
		// ViewX = ObserverOrientation * BigX;
		SmallX[0] = BigX[0];   // Local coordinate system to be displayed by rose of winds
		SmallX[1] = BigX[1];
		SmallX[2] = BigX[2];

		Xlen = K3ToUnity(SmallX, 3);
		if (Xlen > 0.01) {  // show it!
			double SmallY[3], SmallZ[3];
			K3AddUnitProngs(SmallX, SmallY, SmallZ);

			for (int i99 = 0; i99 < 3; i99++) {
				SmallX[i99] *= K3VersorScale;
				SmallY[i99] *= K3VersorScale;
				SmallZ[i99] *= K3VersorScale;
			}

			int VertZeroIs = /* this->*/ vertices().size();


			for (int i = 0; i < 10; i++) {  // looping over the 10 points of arrow;
				tip[0] = RoseStem[0];
				tip[1] = RoseStem[1];
				tip[2] = RoseStem[2];
				for (int i88 = 0; i88 < 3; i88++) {   // looping over three real dimensions
					tip[i88] += K3UnitArrow[0][i] * SmallX[i88] + K3UnitArrow[1][i] * SmallY[i88] + K3UnitArrow[2][i] * SmallZ[i88];
				};
				// i-th point of j-th arrow -if any - ready for insertion
				this->addVertex(CVertex(tip[0], tip[1], tip[2])); // , Kolor);

			}
			// Stem:
			// faces().size;

			// QQ_Ines Spróbować PCA całości i interpretować wyniki

			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 1, VertZeroIs + 2));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 2, VertZeroIs + 3));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 3, VertZeroIs + 4));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 4, VertZeroIs + 1));
			// Base of head:
			this->faces().push_back(CFace(VertZeroIs + 8, VertZeroIs + 5, VertZeroIs + 6));
			this->faces().push_back(CFace(VertZeroIs + 8, VertZeroIs + 6, VertZeroIs + 7));
			// Head:
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 8, VertZeroIs + 7));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 7, VertZeroIs + 6));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 6, VertZeroIs + 5));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 5, VertZeroIs + 8));

			for (int i = 0; i < 10; i++) {  // looping over the 10 faces of arrow;

				this->fcolors().push_back(colorlist[j]);
			};

			// K3MyModel->addChild((CMesh*)this);

		}
		// } if Xlen> ...
	} // } for j=...

	// ThisMesh->addVertex(CVertex(0.0, 0.0, 0.0)); // , Kolor);  AP::WORKSPACE::
}



K3RoseOfWinds::K3RoseOfWinds(MatrixXd K3ViewMat, double* RoseStem, CRGBA* colorlist, double K3VerScale)
// This version for 4D only
// Generate a mesh with the four versors projected into 3D
{
	// int j, j88, i;

	//dp-04-02// K3ListMatrix(DATA_PATH("MojaRurza.txt"), K3ViewMat, "K3ViewMat in rose:");
	// CMesh* this = new CMesh();

#define K3ArrowHeadShare 0.1
#define K3ArrowHeadWidth 0.15
#define K3ArrowStemWidth 0.1
	double a = 1 - K3ArrowHeadShare; // the X at neck
	double b = K3ArrowStemWidth;
	double c = K3ArrowHeadWidth;
	double K3UnitArrow[3][10] = {
	0, a, a, a, a, a, a, a, a, 1,
	0, b, b, -b, -b, c, c, -c, -c, 0,
	0, -b, b, b, -b, -c, c, c, -c, 0 };

	// #define K3VersorScale 200.0  // QQ8 było 200.0. Zmiana na 400 nic dobrego nie zrobila.

	double Xlen2; // to check if vector visible


	double /* BigX[K3Dimensionality], ViewX[K3Dimensionality],*/ LIPA[1000], SmallY[3], SmallX1[3], SmallZ[3], tip[3];


	for (int j = 0; j < 4; j++) {
		// for (int j88 = 0; j88 < K3Dimensionality; j88++) {
		//	BigX[j88] = 0.0;
		// };
		// BigX[j] = 1.0;
		// ViewX = ObserverOrientation * BigX;
		// SmallX[0] = BigX[0];   // Local coordinate system to be displayed by rose of winds
		// SmallX[1] = BigX[1];
		// SmallX[2] = BigX[2];
		SmallX1[0] = K3ViewMat(0, j);   // Local coordinate system to be displayed by rose of winds
		SmallX1[1] = K3ViewMat(1, j);
		SmallX1[2] = K3ViewMat(2, j);

		Xlen2 = SmallX1[0] * SmallX1[0] + SmallX1[1] * SmallX1[1] + SmallX1[2] * SmallX1[2];
		if (Xlen2 > 5.0) {  // show it!
			// double SmallY[3], SmallZ[3];
			K3AddUnitProngs(SmallX1, SmallY, SmallZ);

			for (int i99 = 0; i99 < 3; i99++) {


				SmallY[i99] = 0.077 + SmallY[i99] * K3VerScale;
				SmallX1[i99] = 0.077;// +SmallX1[i99] * K3VerScale; // (K3VerScale * 220.0);
				SmallZ[i99] = 0.077 + SmallZ[i99] * K3VerScale;
			}

			int VertZeroIs = /* this->*/ vertices().size();


			for (int i = 0; i < 10; i++) {  // looping over the 10 points of arrow;
				tip[0] = RoseStem[0];
				tip[1] = RoseStem[1];
				tip[2] = RoseStem[2];
				for (int i88 = 0; i88 < 3; i88++) {   // looping over three real dimensions
					tip[i88] += K3UnitArrow[i88][i] * SmallX1[i88] + K3UnitArrow[1][i] * SmallY[i88] + K3UnitArrow[2][i] * SmallZ[i88];
				};
				// i-th point of j-th arrow -if any - ready for insertion
				this->addVertex(CVertex(tip[0], tip[1], tip[2])); // , Kolor);

			}
			// Stem:
			// faces().size;

			// QQ_Ines Spróbować PCA całości i interpretować wyniki

			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 1, VertZeroIs + 2));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 2, VertZeroIs + 3));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 3, VertZeroIs + 4));
			this->faces().push_back(CFace(VertZeroIs, VertZeroIs + 4, VertZeroIs + 1));
			// Base of head:
			this->faces().push_back(CFace(VertZeroIs + 8, VertZeroIs + 5, VertZeroIs + 6));
			this->faces().push_back(CFace(VertZeroIs + 8, VertZeroIs + 6, VertZeroIs + 7));
			// Head:
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 8, VertZeroIs + 7));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 7, VertZeroIs + 6));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 6, VertZeroIs + 5));
			this->faces().push_back(CFace(VertZeroIs + 9, VertZeroIs + 5, VertZeroIs + 8));

			for (int i = 0; i < 10; i++) {  // looping over the 10 faces of arrow;

				this->fcolors().push_back(colorlist[j]);
			};

			// K3MyModel->addChild((CMesh*)this);

		}
		// } if Xlen> ...
	} // } for j=...

	// ThisMesh->addVertex(CVertex(0.0, 0.0, 0.0)); // , Kolor);  AP::WORKSPACE::
}


// :koniec dodatku LL
