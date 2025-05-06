#include "Avatar3D.h"

#include "Global.h"

#include <Eigen/Dense>
#include <QOpenGl.h>


std::pair<std::vector<Eigen::Vector3f>, std::vector<Eigen::Vector3i>> create_box(float w, float h, float d)
{
    float hw = w / 2;
    float hh = h / 2;
    float hd = d / 2;

    std::vector<Eigen::Vector3f> v = std::vector<Eigen::Vector3f>({
        {-hw, -hh, -hd} ,{hw, -hh, -hd}, {hw, hh, -hd}, {-hw, hh, -hd},
        {-hw, -hh, hd}, {hw, -hh, hd}, {hw, hh, hd}, {-hw, hh, hd} });

    std::vector<Eigen::Vector3i> f = std::vector<Eigen::Vector3i>({
        {0, 1, 2}, {2, 3, 0},
        {4, 5, 6}, {6, 7, 4},
        {0, 4, 7}, {7, 3, 0},
        {1, 5, 6}, {6, 2, 1},
        {3, 2, 6}, {6, 7, 3},
        {0, 1, 5}, {5, 4, 0} });

    return { v, f };
}

Eigen::Vector3f skin_color_from_feature(float val)
{
    float value = val;

    if (value < -1.0) value = -1.0;
    else if (value > 1.0) value = 1.0;

    Eigen::Vector3f light_skin(1.0, 0.85, 0.7);
    Eigen::Vector3f dark_skin(0.5, 0.35, 0.2);
    float t = (value + 1.0) / 2.0;  // [-1, 1] ->[0, 1]

    return ((1 - t) * dark_skin + t * light_skin);
}


void AvatarPart::draw(Eigen::Vector3f offset)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

    glColor3f(m_color[0], m_color[1], m_color[2]);
    glBegin(GL_TRIANGLES);

    for (Eigen::Vector3i face : m_faces)
        for (int j = 0; j < 3; j++)
        {
            int idx = face[j];
            auto v = m_vertices[idx] + m_position + offset;
            glVertex3f(v[0], v[1], v[2]);
        }
    glEnd();

    glDisable(GL_COLOR_MATERIAL);

    glPopAttrib();
}
