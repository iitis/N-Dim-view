#include "K3Totem.h"

#include "K3Helpers.h"

//void add_cube(CMesh* ThisMesh, CRGBA Kolor, const double scale[] = {1.0,1.0,1.0})
//{
//	int vof = ThisMesh->vertices().size();
//
//	std::vector<CVertex> vcs({
//		CVertex(0.0, 0.0, 0.0), CVertex(0.0, 0.0, 1.0), CVertex(0.0, 1.0, 0.0), CVertex(0.0, 1.0, 1.0),
//		CVertex(1.0, 0.0, 0.0), CVertex(1.0, 0.0, 1.0), CVertex(1.0, 1.0, 0.0), CVertex(1.0, 1.0, 1.0) });
//	
//	ThisMesh->vertices().insert(ThisMesh->vertices().end(), vcs.begin(), vcs.end());
//
//	std::vector<CFace> fcs({
//		CFace(1, 5, 3),CFace(7, 3, 5),
//		CFace(4, 6, 5),CFace(7, 5, 6),
//		CFace(0, 2, 4),CFace(6, 4, 2),
//		CFace(1, 3, 0),CFace(2, 0, 3),
//		CFace(2, 3, 6),CFace(7, 6, 3),
//		CFace(0, 4, 1),CFace(5, 1, 4) });
//
//	for (auto& f : fcs) {
//		f[0] += vof;
//		f[1] += vof;
//		f[2] += vof;
//	}
//
//	ThisMesh->faces().insert(ThisMesh->faces().end(), fcs.begin(), fcs.end());
//
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
//}

void K3Totem::K3FillMeshToUnitCube(CMesh* ThisMesh, CRGBA Kolor) { // double l, double a) {
	//CMesh* K3UnitCube = new CMesh();
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


std::shared_ptr<K3ChernoffFace> K3Totem::make_face(std::vector<std::optional<double>> values, CMesh* korpus)
{
	//int nnn = K3HyperLook.rows();
	//double k3xx[] = { 0.4,0.3,0.8,0.9,0.7,0.6,0.5,0.4,0.7, 0.9 };
	//for (int i = 0; i < std::min(10, nnn); i++) { // size(K3HyperLook,1)
	//	k3xx[i] = K3HyperLook[i];
	//}

	double cx, cy;
	//K3ChernoffFace* K3Cher3;   // 9781 1159
	std::shared_ptr<K3ChernoffFace> K3Cher3 = std::make_shared<K3ChernoffFace>(cx = 600, cy = 800, values);
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

	K3Cher3->setLabel("face");

	return K3Cher3;
}


K3Totem::K3Totem(Eigen::VectorXd K3HyperSpot, Eigen::VectorXd K3HyperLook) {
	// Create a totem whose position and appearance represents data.
	// Assume the given K3HyperSpot = DataPoint*Observer ,
	// i.e. has been projected into the observation space.
	std::shared_ptr<CMesh> korpus = std::make_shared<CMesh>();
	std::shared_ptr<CMesh> LeftArm = std::make_shared<CMesh>();

	// UI::MESSAGEBOX::error(L"Totem potem");

	//Eigen::Matrix4d K3Position;
	//K3Position <<
	//	0.51, 0.0, 0.0, K3HyperSpot[0],
	//	0.0, 0.5, 0.0, K3HyperSpot[1],
	//	0.0, 0.0, 0.5, K3HyperSpot[2],
	//	0.0, 0.0, 0.0, 0.2;

	Eigen::Matrix4d K3Position = Eigen::Matrix4d::Zero();
	K3Position.diagonal() = Eigen::Vector4d(0.51, 0.5, 0.5, 0.2);
	K3Position.topRightCorner(3, 1) = K3HyperSpot.topRows(3);

	K3FillMeshToUnitCube(korpus.get(), CRGBA(0.2f, 0.20f, 1.0f, 1.0f)); // K3_color(K3HyperLook[19] / 1.0, 1.0)); // CRGBA(1.0f, 1.0f, 0.0f, 0.6f));

	// Body color
	// CAnnotationPoint *K3Navel = new CAnnotationPoint(0.5, .5, 1.0);

	Eigen::Matrix4d K3BodyShape = Eigen::Matrix4d::Zero();
	K3BodyShape.diagonal() = Eigen::Vector4d(1.0, 2.0, 0.6, 1.0);

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

	korpus->setLabel("body");
	this->addChild(korpus);
	// this->addChild(K3Navel);
	// this->addChild(K3Ad);

	// Make left arm:
		// QQSize size of the arm determined here:

	K3FillMeshToUnitCube(LeftArm.get(), CRGBA(0.0f, 1.0f, 0.0f, 1.0f));

	Eigen::Matrix4d K3LeftArmPosition;

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


	// this->addChild(LeftArm);


	// delete korpus;?
	// delete LeftArm;

	// Now let's put a face on it:

	// tymczasowo:
	std::vector<std::optional<double>> values( { 0.4, 0.3, 0.8, 0.9, 0.7, 0.6, 0.5, 0.4, 0.7, 0.9 } );
	// std::vector<std::optional<double>> values(10);

	for (int i = 0; i < std::min(10, static_cast<int>(K3HyperLook.rows())); i++)
		values[i] = K3HyperLook[i];

	this->addChild( make_face(values, korpus.get() ) );

};
// end of K3Totem
