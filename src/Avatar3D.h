#pragma once
#include <Object.h>
#include <Eigen/Dense>
#include "dll_global.h"

class DPVISION_DLL_API AvatarPart
{
public:
    std::vector<Eigen::Vector3f> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;
    Eigen::Vector3f m_color;
    Eigen::Vector3f m_position;

    AvatarPart();
    void draw(Eigen::Vector3f offset);
    void set(std::vector<Eigen::Vector3f> vertices, std::vector<Eigen::Vector3i> faces, Eigen::Vector3f color);
    void translate(float dx, float dy, float dz);
};



class Head : public AvatarPart
{
public:
    Head(double feature_0, double feature_6);
};


class Nose : public AvatarPart
{
public:
    Nose(double feature_1);
};

class Eye : public AvatarPart
{
public:
    Eye(double feature_5, bool left = true);
};


class Mouth : public AvatarPart
{
public:
    Mouth(double feature_2, double feature_3);
    void create_curved_mouth(double width, double curve, double thickness, int segments,
        std::vector<Eigen::Vector3f>& vertices, std::vector<Eigen::Vector3i>& faces);
};


class Hair : public AvatarPart
{
public:
    Hair(double feature_4);
};


class DPVISION_DLL_API Avatar3D : public CObject
{
    std::vector<std::shared_ptr<AvatarPart>> m_parts;

    std::vector<Eigen::Vector3f> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;
    std::vector<Eigen::Vector3f> m_colors;

    Eigen::Vector3f m_position;

public:
    Avatar3D(const Eigen::VectorXd& features, const Eigen::Vector3d& pos = Eigen::Vector3d::Zero(4));

    //virtual inline int type() { return CBaseObject::Type::NEWTYPES + 7; };

    void translate(const Eigen::Vector3f& delta);
    void buildMesh();
    void draw() const;
    virtual void renderSelf() override;
    void buildCombinedMesh();
};


class QMouseEvent;
class QWheelEvent;

class PropWidget;

class DPVISION_DLL_API SwarmOfAvatars3D : public CObject
{
public:
    using CObject::CObject;
    virtual bool mousePressEvent(QMouseEvent* event) override;
    virtual bool mouseReleaseEvent(QMouseEvent* event) override;
    virtual bool mouseMoveEvent(QMouseEvent* event) override;
    virtual bool wheelEvent(QWheelEvent* event) override;

    //inline virtual bool has_prop_widget_defined() override { return true; };

    //QVector<PropWidget*> create_propwidget_and_get_subwidgets();

};

//#include "propWidget.h"
//#include <QVBoxLayout>
//
//class DPVISION_DLL_API PropSwarmOfAvatars3D : public PropWidget
//{
//    Q_OBJECT
//
//    SwarmOfAvatars3D* obj;
//
//public:
//    explicit PropSwarmOfAvatars3D(SwarmOfAvatars3D* swarm, QWidget* parent = 0) : PropWidget(parent) {
//        obj = swarm;
//
//        treeItemLabel = "SwarmOfAvatars3D properties";
//
//        //ui.setupUi(this);
//
//        //connect(ui.cloud, &QGroupBox::toggled, this, &PropPointCloud::adjustGroupHeight);
//    };
//
//
//    ~PropSwarmOfAvatars3D() {};
//
//
//    static PropWidget* create(SwarmOfAvatars3D* m, QWidget* parent = 0) {
//        PropWidget* widget = new PropWidget(parent);
//        QVBoxLayout* layout = new QVBoxLayout(widget);
//
//        //layout->addWidget(new PropBaseObject(m, widget));
//        layout->addWidget(new PropSwarmOfAvatars3D(m, widget));
//
//        widget->resize(layout->sizeHint());
//        widget->setMinimumSize(layout->sizeHint());
//        widget->setMaximumSize(layout->sizeHint());
//
//        widget->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
//
//        return widget;
//    };
//    
//    static QVector<PropWidget*> create_and_get_subwidgets(CBaseObject* obj) {
//        SwarmOfAvatars3D* m = (SwarmOfAvatars3D*)obj;
//        return QVector<PropWidget*>({
//            //new PropBaseObject(m),
//            new PropSwarmOfAvatars3D(m)
//        });
//
//        return QVector<PropWidget*>();
//    };
//
//    virtual void updateProperties() override {};
//
//public slots:
//
//private:
//
//};
