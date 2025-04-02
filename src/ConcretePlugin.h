#pragma once

#include "dll_global.h"
#include "PluginInterface.h"

#include <Eigen/Geometry>

#include "CsvColumnAssignmentDialog.h"
#include "CsvReader.h"

class CModel3D;
class QSettings;

#define PNG_DIR "d:/K3/Wielowymiar/"

class DPVISION_DLL_API ConcretePlugin : public QObject, public PluginInterface
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "dpVision.PluginInterface" FILE "metadata.json")
	Q_INTERFACES(PluginInterface)

	CsvReader::CsvData current_data;

	Eigen::MatrixXd *K3_IObs, K3BoiledData;

	QVector<ColumnAssignment> current_assignment;

	// Kazda linia reprezentuje grupê, pierwszy string to jej nazwa
	// nastêpnie w {} s¹ etykiety zdefiniwanych przez ciebie wymiarów w tej grupie
	// Jesli {} jest puste tzn, ¿e grupa nie moze mieæ etykiet (np.: nienazwane)
	// UWAGA - to nie sa nazwy wymiarów zwi¹zanych z plikiem Ÿród³owym (np. wino)
	// tylko to jak je bêdziesz u¿ywa³ podczas wizualizowalizacji
	// powi¹zanie tego z np. winem dzieje sie dopiero w okienku dialogowym
	// dlatego jest to niezale¿ne od rodzaju danych
	// Mo¿esz sobie to dowolnie zmodyfikowaæ, np dodac nowe elementy
	QVector<GroupDefinition> groupDefs = {
		{ "unnamed", {} },
		{ "spacial", { "X", "Y", "Z", "T" } },
		{ "visual", { "Skin_C", "Hair_C", "Eye_S", "Nose_L", "Mouth_W", "Smile", "Frown", "Hair_L", "Face_Elong", "Iris_C"}},
	};

	QString png_save_dir;


	class Registry {
	public:
		static void storeAssignments(QSettings* settings, const QString& filePath, const QVector<QPair<int, int>>& vec);
		static QVector<QPair<int, int>> loadAssignments(QSettings* settings, const QString& filePath);
		static void cleanupOldAssignments(QSettings* settings, int daysOld);
	};

	Eigen::MatrixXd K3_Get_PCA_Funnel(Eigen::MatrixXd X, int nd);
	void K3AddMyCloud(CModel3D* K3MyModel, Eigen::MatrixXd K3ObsCloud, Eigen::MatrixXd K3ViewMat, double K3Toler);
	int K3FormProjectionMatrix(const Eigen::MatrixXd& RawData);
	void setDatasetLabel();
	void setAssignmentLabel();
	void delete_old_screenshots(const QString& pattern);
	void AssignGroups();
	void K3CzteroPajak(double h05);
	void K3Krata(int Nkrat, int Mkrat);
	void K3Display();
	void K3Dance(double grand_scale);
	void K3LoadDataset();

public:
    ConcretePlugin(void);
	~ConcretePlugin(void) override;

	/* ============= EVENT HANDLING BEGIN =========================*/
	virtual void onLoad() override;
	virtual void onButton(std::wstring name) override;
	/* ============= EVENT HANDLING END =========================*/
};
