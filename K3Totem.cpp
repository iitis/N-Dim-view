#include "K3Totem.h"
#include "K3ChernoffFace.h"

#include "K3Helpers.h"

K3Totem::K3Totem(Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook) {
	// Create a totem whose position and appearance represents data.
	// Assume the given K3HyperSpot = DataPoint*Observer ,
	// i.e. has been projected into the observation space.
	CMesh* korpus = new CMesh;
	CMesh* LeftArm = new CMesh;

	// UI::MESSAGEBOX::error(L"Totem potem");

	Eigen::Matrix4d K3BodyShape, K3Position, K3LeftArmPosition;
	K3Position(0, 0) = .51; //  1.0;   // QQ7 SCALING
	K3Position(0, 1) = 0.0;
	K3Position(0, 2) = 0.0;
	K3Position(0, 3) = K3HyperSpot[0];

	K3Position(1, 0) = 0.0;
	K3Position(1, 1) = .5; // 1.0;
	K3Position(1, 2) = 0.0;
	K3Position(1, 3) = K3HyperSpot[1];

	K3Position(2, 0) = 0.0;
	K3Position(2, 1) = 0.0;
	K3Position(2, 2) = .5; // 1.0;
	K3Position(2, 3) = K3HyperSpot[2];

	K3Position(3, 0) = 0.0;
	K3Position(3, 1) = 0.0;
	K3Position(3, 2) = 0.0;
	K3Position(3, 3) = 0.2;    // QQVISUS  (was 1.0)

	\
		/* Od DP:
		* akurat w wymienionym przypadku mo�na jeszcze pro�ciej:

	Eigen::Matrix4d K3Position = Eigen::Matrix4d::Identity();

	albo:
	Eigen::Matrix4d K3Position;
	K3Position.setIdentity();

	ale nie zawsze mamy takie oczywiste macierze :)
		*/

		K3FillMeshToUnitCube(korpus, K3_color(K3HyperLook[19] / 1.0, 1.0)); // CRGBA(1.0f, 1.0f, 0.0f, 0.6f));
	// Body color
	// CAnnotationPoint *K3Navel = new CAnnotationPoint(0.5, .5, 1.0);


	K3BodyShape.setIdentity();
	K3BodyShape(0, 0) = 1.0;
	K3BodyShape(1, 1) = 2.0;
	K3BodyShape(2, 2) = 0.6;
	for (int i = 0; i < korpus->vertices().size(); i++)
	{
		// QQSize bodyshape in the line below may impact size:
		korpus->vertices()[i] = korpus->vertices()[i].transformByMatrix(K3BodyShape);
	}
	// CAnnotation* K3Ad = new CAnnotation;
	// K3Ad->setLabel("LAx");
	// K3Navel->addAnnotation(K3Ad);

	CPoint3<float> K3LeftShoulder; // = [0.7, 0.7, 0];
	K3LeftShoulder[0] = 0.7;
	K3LeftShoulder[1] = 1.5;
	K3LeftShoulder[2] = 0.05;

	for (int i = 0; i < korpus->vertices().size(); i++)
	{
		korpus->vertices()[i] = korpus->vertices()[i].transformByMatrix(K3Position);
	}
	K3LeftShoulder.transformByMatrix(K3Position);


	this->addChild(korpus);
	// this->addChild(K3Navel);
	// this->addChild(K3Ad);

	// Make left arm:
		// QQSize size of the arm determined here:

	K3FillMeshToUnitCube(LeftArm, CRGBA(0.0f, 1.0f, 0.0f, 1.0f));
	K3LeftArmPosition.setIdentity(); // = Eigen::Matrix4d::Identity();
	K3LeftArmPosition(0, 0) = 1.7; // arm is long  (was 1.0 to 2.3)
	K3LeftArmPosition(1, 1) = .2; // and slim
	K3LeftArmPosition(2, 2) = .2; // and slim
	K3LeftArmPosition(0, 3) = -0.5; // to get pivot in O, shift in x ...
	K3LeftArmPosition(1, 3) = -0.1; // ...and y...
	K3LeftArmPosition(2, 3) = 0.1; // ...and z.
	for (int i = 0; i < LeftArm->vertices().size(); i++)
	{
		LeftArm->vertices()[i] = LeftArm->vertices()[i].transformByMatrix(K3LeftArmPosition);
	}


	// Now raise arm according to parameter:
	double K3Alpha = K3HyperLook[10]; // qq? brow
	K3LeftArmPosition.setIdentity(); // = Eigen::Matrix4d::Identity();
	K3LeftArmPosition(0, 0) = cos(K3Alpha);
	K3LeftArmPosition(0, 1) = -sin(K3Alpha);
	K3LeftArmPosition(1, 0) = sin(K3Alpha);
	K3LeftArmPosition(1, 1) = cos(K3Alpha);
	// Now shift to shoulder:
	K3LeftArmPosition(0, 3) += K3LeftShoulder[0];
	K3LeftArmPosition(1, 3) += K3LeftShoulder[1];
	K3LeftArmPosition(2, 3) += ((double)(K3LeftShoulder[2]));

	for (int i = 0; i < LeftArm->vertices().size(); i++)
	{
		LeftArm->vertices()[i] = LeftArm->vertices()[i].transformByMatrix(K3LeftArmPosition);
	}

	// Now move arm to trunk:
	K3LeftArmPosition = K3Position;
	for (int i = 0; i < LeftArm->vertices().size(); i++)
	{
		LeftArm->vertices()[i] = LeftArm->vertices()[i].transformByMatrix(K3LeftArmPosition);
	}


	this->addChild(LeftArm);


	// delete korpus;?
	// delete LeftArm;

	// Now let's put a face on it:
	if (1) {   // NO FACE!
		{
			int nnn = K3HyperLook.rows();
			double k3xx[] = { 0.4,0.3,0.8,0.9,0.7,0.6,0.5,0.4,0.7, 0.9 };
			for (int i = 0; i < std::min(10, nnn); i++) { // size(K3HyperLook,1)
				k3xx[i] = K3HyperLook[i];
			}

			double cx, cy;
			//K3ChernoffFace* K3Cher3;   // 9781 1159
			K3ChernoffFace* K3Cher3 = new K3ChernoffFace(cx = 600, cy = 800, k3xx);
			//CImage* CimCher3 = (CImage*)K3Cher3;
			//CTransform K3ImagePos;
			//K3ImagePos = CimCher3->getTransform();


			double K3D[16] = { 2.0,0.0,0.0,0.0, 0.0,1.0,0.0,0.0,  0.0,0.0,1.0,44.0,  0.0,0.0,10.0,1.0 };

			// korpus->vertices()[i]

			Eigen::Matrix4Xd K3SourceFaceFrame(4, 5), K3DestFaceFrame(4, 5);


			Eigen::Matrix4d K3FaceToHead;
			K3SourceFaceFrame(0, 0) = cx / 40.;
			K3SourceFaceFrame(0, 1) = -cx / 40.;
			K3SourceFaceFrame(0, 2) = -cx / 40.;
			// K3SourceFaceFrame(0, 3) = cx / 40.;
			K3SourceFaceFrame(0, 3) = 0.0;

			K3SourceFaceFrame(1, 0) = cy / 40.;
			K3SourceFaceFrame(1, 1) = cy / 40.;
			K3SourceFaceFrame(1, 2) = -cy / 40.;
			// K3SourceFaceFrame(1, 3) = -cy / 40.;
			K3SourceFaceFrame(1, 3) = 0.0;

			K3SourceFaceFrame(2, 0) = 0.0;
			K3SourceFaceFrame(2, 1) = 0.0;
			K3SourceFaceFrame(2, 2) = 0.0;
			// K3SourceFaceFrame(2, 3) = 0.0;
			K3SourceFaceFrame(2, 3) = 10.0;

			K3SourceFaceFrame(3, 0) = 1.0;
			K3SourceFaceFrame(3, 1) = 1.0;
			K3SourceFaceFrame(3, 2) = 1.0;
			// K3SourceFaceFrame(3, 3) = 1.0;
			K3SourceFaceFrame(3, 3) = 1.0;

			double drawbar[3];
			double L1[3], L2[3];
			for (int i88 = 0; i88 < 3; i88++) {
				L1[i88] = korpus->vertices()[5][i88] / korpus->vertices()[5][3] - korpus->vertices()[1][i88] / korpus->vertices()[1][3];
				L2[i88] = korpus->vertices()[3][i88] / korpus->vertices()[3][3] - korpus->vertices()[1][i88] / korpus->vertices()[1][3];
			}
			drawbar[0] = L1[1] * L2[2] - L2[1] * L1[2];
			drawbar[1] = L1[2] * L2[0] - L2[2] * L1[0];
			drawbar[2] = L1[0] * L2[1] - L2[0] * L1[1];


			K3DestFaceFrame(0, 0) = korpus->vertices()[7][0];
			K3DestFaceFrame(0, 1) = korpus->vertices()[3][0];
			K3DestFaceFrame(0, 2) = korpus->vertices()[1][0];
			// K3DestFaceFrame(0, 3) = korpus->vertices()[5][0];
			K3DestFaceFrame(0, 3) = (korpus->vertices()[0][0] + korpus->vertices()[2][0] + korpus->vertices()[4][0] + korpus->vertices()[6][0]) / 4.0; // +drawbar[0];

			K3DestFaceFrame(1, 0) = korpus->vertices()[7][1];
			K3DestFaceFrame(1, 1) = korpus->vertices()[3][1];
			K3DestFaceFrame(1, 2) = korpus->vertices()[1][1];
			// K3DestFaceFrame(1, 3) = korpus->vertices()[5][1];
			K3DestFaceFrame(1, 3) = (korpus->vertices()[0][1] + korpus->vertices()[2][1] + korpus->vertices()[4][1] + korpus->vertices()[6][1]) / 4.0; // + drawbar[1];

			K3DestFaceFrame(2, 0) = .01 + korpus->vertices()[7][2];
			K3DestFaceFrame(2, 1) = .01 + korpus->vertices()[3][2];
			K3DestFaceFrame(2, 2) = .01 + korpus->vertices()[1][2];
			// K3DestFaceFrame(2, 3) = korpus->vertices()[5][2];
			K3DestFaceFrame(2, 3) = .01 + (korpus->vertices()[0][2] + korpus->vertices()[2][2] + korpus->vertices()[4][2] + korpus->vertices()[6][2]) / 4.0; // + drawbar[2];

			K3DestFaceFrame(3, 0) = 1.0;
			K3DestFaceFrame(3, 1) = 1.0;
			K3DestFaceFrame(3, 2) = 1.0;
			// K3DestFaceFrame(3, 3) = 1.0;
			K3DestFaceFrame(3, 3) = 1.0;

			K3FindTransform(K3SourceFaceFrame, K3DestFaceFrame, &K3FaceToHead, 4); // ( K3DestFaceFrame.inv   .inverse());
			double K3FaceToHeadRaw[16];
			for (int j88 = 0; j88 < 4; j88++) {
				for (int i88 = 0; i88 < 4; i88++) {
					K3FaceToHeadRaw[4 * j88 + i88] = K3FaceToHead(j88, i88);   // QQ: Trying transposition!
				}

			}

			// K3ImagePos.fromRowMatrixD(K3FaceToHeadRaw);
			// K3Cher3->getTransform().fromRowMatrixD(K3D);
			K3Cher3->getTransform().fromRowMatrixD(K3FaceToHeadRaw); // K3FaceToHeadRaw or K3D

			this->addChild((CModel3D*)K3Cher3); // QQ Wlknoc blok
		}
	}
};
// end of K3Totem
