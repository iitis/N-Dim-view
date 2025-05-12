#include "SwarmOfAvatars3D.h"

#include "Global.h"

#include <Eigen/Dense>
#include <QOpenGl.h>

bool SwarmOfAvatars3D::mousePressEvent(QMouseEvent* event)
{
    qWarning() << "(Swarm) mouse pressed";
    return true;
}

bool SwarmOfAvatars3D::mouseReleaseEvent(QMouseEvent* event)
{
    qWarning() << "(Swarm) mouse released";
    return true;
}

bool SwarmOfAvatars3D::mouseMoveEvent(QMouseEvent* event)
{
    qWarning() << "(Swarm) mouse moved";
    return true;
}

bool SwarmOfAvatars3D::wheelEvent(QWheelEvent* event)
{
    qWarning() << "(Swarm) mouse wheel event";
    return true;
}

#include <QtWidgets>

QWidget* SwarmOfAvatars3D::prop_widget() {
    m_prop_widget = new SwarmWidget();

    return (QWidget*)m_prop_widget;
}

SwarmWidget::SwarmWidget(QWidget* parent) : QWidget(parent)
{
    QVBoxLayout* layout = new QVBoxLayout((QWidget *)this);

    layout->addWidget(new QLabel(QString::fromUtf8("MIEJSCE NA TWOJĄ REKLAMĘ")));

    this->resize(layout->sizeHint());
    this->setMinimumSize(layout->sizeHint());
    this->setMaximumSize(layout->sizeHint());

    this->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);

}
