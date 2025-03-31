#pragma once

#include "dll_global.h"

#include "PluginInterface.h"
#include "Mesh.h"
#include "Image.h"
#include <Eigen/Geometry>

#include "CsvColumnAssignmentDialog.h"

class K3Totem;

class DPVISION_DLL_API ConcretePlugin : public QObject, public PluginInterface
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "dpVision.PluginInterface" FILE "metadata.json")
	Q_INTERFACES(PluginInterface)

	bool m_picking;

	Eigen::MatrixXd current_data_matrix, *K3_IObs, K3BoiledData;
	QVector<ColumnAssignment> current_assignment;



public:
    ConcretePlugin(void);
	~ConcretePlugin(void) override;

	void K3Display(CModel3D* o, float d);
	// Added by LL:
	Eigen::MatrixXd K3_Get_PCA_Funnel(Eigen::MatrixXd X, int nd);
	int K3FormProjectionMatrix(Eigen::MatrixXd* RawData);

	// :added by LL


	/* ============= EVENT HANDLING BEGIN =========================*/

	// called on doubleclick on plugin name
	//virtual void onRun(void) override {};

	// called when plugin is loaded
	virtual void onLoad() override;

	// The way to read in data from csv file:
	// https://java2blog.com/read-csv-file-in-cpp/#How_To_Read_A_CSV_File_In_C



	// called before plugin unload
	//virtual void onUnload() override {};

	// called when you select plugin
	//virtual void onActivate() override {};

	// called when you select another plugin
	//virtual void onDeactivate() override {};

	// called if button on PluginPanel is pressed
	virtual void onButton(std::wstring name) override;

	// called if plugin is active and ctrl+MouseClick is done on any model
	//virtual bool onMousePick(Plugin::PickEvent pickEvent) override;

	/* ============= EVENT HANDLING END =========================*/

};
