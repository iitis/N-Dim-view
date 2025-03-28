#include "K3Arrow.h"

#include "K3Helpers.h"

#include "RGBA.h"

K3Arrow::K3Arrow(double A[3], double B[3], double R, CRGBA* colour) : CMesh(nullptr) {
	// Make 3D arrow from A to B, with thickness R and colour colour

	double D[3], L[3], T[3];  // D = free vector; L,T = fletchings (random orientation)
	int i;
#define K3ArrowHeadShare 0.35
	// used to be 0.35
#define K3ArrowHeadWidth 0.15
#define K3ArrowStemWidth 0.05
	double a = 1 - K3ArrowHeadShare; // the X at neck
	double b = K3ArrowStemWidth;
	double c = K3ArrowHeadWidth;
	double K3UnitArrow[3][10] = {
	0, a, a, a, a, a, a, a, a, 1,
	0, b, b, -b, -b, c, c, -c, -c, 0,
	0, -b, b, b, -b, -c, c, c, -c, 0 };
	double SmallX[3], SmallY[3], SmallZ[3], tip[3];

	for (i = 0; i < 3; i++) {
		D[i] = B[i] - A[i];
	};
	K3AddUnitProngs(D, L, T);

	SmallX[0] = D[0];
	SmallX[1] = D[1];
	SmallX[2] = D[2];

	// Xlen = K3ToUnity(SmallX, 3);
	if (2.0 > 0.01) {  // show it!
		double SmallY[3], SmallZ[3];
		K3AddUnitProngs(SmallX, SmallY, SmallZ);
		for (int i88 = 0; i88 < 3; i88++) {
			SmallY[i88] *= R;
			SmallZ[i88] *= R;
		}

		int VertZeroIs = this->vertices().size();

		for (int i = 0; i < 10; i++) {  // looping over the 10 points of arrow;
			tip[0] = A[0];
			tip[1] = A[1];
			tip[2] = A[2];
			for (int i88 = 0; i88 < 3; i88++) {   // looping over three real dimensions
				tip[i88] += K3UnitArrow[0][i] * SmallX[i88] + K3UnitArrow[1][i] * SmallY[i88] + K3UnitArrow[2][i] * SmallZ[i88];
			};
			// i-th point of j-th arrow -if any - ready for insertion
			this->addVertex(CVertex(tip[0], tip[1], tip[2])); // , Kolor);

		}
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

			this->fcolors().push_back(*colour);
		};

	}
};
