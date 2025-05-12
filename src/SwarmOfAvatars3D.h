#pragma once
#include <Object.h>
#include <Eigen/Dense>
#include "dll_global.h"
#include <QWidget>

class QMouseEvent;
class QWheelEvent;


class SwarmWidget : public QWidget
{
    Q_OBJECT

public:
    explicit SwarmWidget(QWidget* parent = 0);
};

class DPVISION_DLL_API SwarmOfAvatars3D : public CObject
{
    SwarmWidget* m_prop_widget;

public:
    using CObject::CObject;
    virtual bool mousePressEvent(QMouseEvent* event) override;
    virtual bool mouseReleaseEvent(QMouseEvent* event) override;
    virtual bool mouseMoveEvent(QMouseEvent* event) override;
    virtual bool wheelEvent(QWheelEvent* event) override;

    inline virtual bool has_prop_widget() override { return true; };
    virtual QWidget* prop_widget() override;
    inline virtual void prop_widget_update() override {};
};

