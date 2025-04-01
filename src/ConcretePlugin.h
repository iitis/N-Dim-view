#pragma once

#include "dll_global.h"
#include "PluginInterface.h"

#include <Eigen/Geometry>

#include "CsvColumnAssignmentDialog.h"
#include "CsvReader.h"

class CModel3D;

class DPVISION_DLL_API ConcretePlugin : public QObject, public PluginInterface
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "dpVision.PluginInterface" FILE "metadata.json")
	Q_INTERFACES(PluginInterface)

	CsvReader::CsvData current_data;
	//Eigen::MatrixXd current_data_matrix;
	Eigen::MatrixXd *K3_IObs, K3BoiledData;
	QVector<ColumnAssignment> current_assignment;

public:
    ConcretePlugin(void);
	~ConcretePlugin(void) override;

	void K3Display(CModel3D* o, float d);
	Eigen::MatrixXd K3_Get_PCA_Funnel(Eigen::MatrixXd X, int nd);
	int K3FormProjectionMatrix(const Eigen::MatrixXd &RawData);
	//int K3FormProjectionMatrix(Eigen::MatrixXd* RawData);

	void AssignGroups();
	void K3CzteroPajak(double h05);
	void K3Krata(int Nkrat, int Mkrat);
	void K3Display();
	void K3Dance(double grand_scale);


	/* ============= EVENT HANDLING BEGIN =========================*/
	virtual void onLoad() override;
	virtual void onButton(std::wstring name) override;
	/* ============= EVENT HANDLING END =========================*/
};
