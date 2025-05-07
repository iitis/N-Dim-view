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
    //Eigen::Vector3f light_skin(1.0, 1.0, 0.0);
    //Eigen::Vector3f dark_skin(0.0, 1.0, 0.0);
    float t = (value + 1.0) / 2.0;  // [-1, 1] ->[0, 1]

    return ((1 - t) * dark_skin + t * light_skin);
}


inline AvatarPart::AvatarPart()
{
    m_vertices = std::vector<Eigen::Vector3f>();
    m_faces = std::vector<Eigen::Vector3i>();
    m_color = { 0,0,0 };
    m_position = { 0, 0, 0 };
}

inline void AvatarPart::draw(Eigen::Vector3f offset)
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

inline void AvatarPart::set(std::vector<Eigen::Vector3f> vertices, std::vector<Eigen::Vector3i> faces, Eigen::Vector3f color)
{
    m_vertices = vertices;
    m_faces = faces;
    m_color = color;
    m_position = { 0, 0, 0 };
}

inline void AvatarPart::translate(float dx, float dy, float dz)
{
    m_position += Eigen::Vector3f(dx, dy, dz);
}

inline Avatar3D::Avatar3D(const Eigen::VectorXd& features, const Eigen::Vector3d& pos) : CObject()
{
    Eigen::VectorXd feat = features; // (maxAbs.array().abs() > 1e-12).select(features.array() / maxAbs.array(), features.array());

    m_parts.push_back(std::make_shared<Head>(feat[8], feat[0]));
    m_parts.push_back(std::make_shared<Hair>(feat[1]));
    m_parts.push_back(std::make_shared<Nose>(feat[3]));
    m_parts.push_back(std::make_shared<Eye>(feat[2], true));
    m_parts.push_back(std::make_shared<Eye>(feat[2], false));
    m_parts.push_back(std::make_shared<Mouth>(feat[4], feat[5]));

    m_position = pos.cast<float>();
    //buildMesh();
}

inline void Avatar3D::translate(const Eigen::Vector3f& delta) {
    m_position += delta;
}

inline void Avatar3D::buildMesh() {
    buildCombinedMesh();
    // VBO creation omitted for brevity
}

inline void Avatar3D::draw() const {
    for (const auto& part : m_parts) {
        part->draw(m_position);
    }
}

inline void Avatar3D::renderSelf()
{
    draw();
}

inline void Avatar3D::buildCombinedMesh() {
    m_vertices.clear();
    m_faces.clear();
    m_colors.clear();
    int vertexOffset = 0;

    for (const auto& part : m_parts) {
        const auto& v = part->m_vertices;
        const auto& f = part->m_faces;
        const auto& c = part->m_color;
        Eigen::Vector3f partPos = part->m_position;

        for (const auto& vert : v)
            m_vertices.push_back(vert + partPos);
        for (const auto& face : f)
            m_faces.push_back(face + Eigen::Vector3i::Constant(vertexOffset));
        for (size_t i = 0; i < f.size(); ++i)
            m_colors.push_back(c);

        vertexOffset += static_cast<int>(v.size());
    }
}

inline Hair::Hair(double feature_4)
{
    double hair_color = (feature_4 + 1.0) / 2.0;
    Eigen::Vector3f c(hair_color, hair_color, hair_color);
    auto [v, f] = create_box(1.2, 1.0, 0.6);
    set(v, f, c);
    translate(0, 0.4, -0.1);
}

inline Mouth::Mouth(double feature_2, double feature_3)
{
    double width = 0.6 + feature_2 / 2.0;
    double curve = 0.0 + feature_3 / 5.0;
    double thickness = 0.1;
    int segments = 6;

    create_curved_mouth(width, curve, thickness, segments, m_vertices, m_faces);
    m_color = { 1.0, 0.0, 0.0 };

    translate(0, -0.3, 0.31);
}

inline void Mouth::create_curved_mouth(double width, double curve, double thickness, int segments, std::vector<Eigen::Vector3f>& vertices, std::vector<Eigen::Vector3i>& faces)
{

    vertices.clear();
    faces.clear();

    std::vector<float> x_vals(segments + 1);
    std::vector<float> y_vals(segments + 1);

    for (int i = 0; i <= segments; ++i) {
        float t = static_cast<float>(i) / segments;
        x_vals[i] = -width / 2.0 + t * width;
        y_vals[i] = -curve * (1.0 - std::pow((2.0 * x_vals[i] / width), 2));
    }

    int vert_offset = 0;
    float z = 0.0f;

    for (int i = 0; i < segments; ++i) {
        Eigen::Vector3f p1(x_vals[i], y_vals[i], z);
        Eigen::Vector3f p2(x_vals[i + 1], y_vals[i + 1], z);

        Eigen::Vector3f direction = p2 - p1;
        float length = direction.norm();
        if (length < 1e-6) continue;

        direction.normalize();
        Eigen::Vector3f offset(-direction.y(), direction.x(), 0.0);
        offset *= (thickness / 2.0);

        Eigen::Vector3f v0 = p1 + offset;
        Eigen::Vector3f v1 = p1 - offset;
        Eigen::Vector3f v2 = p2 - offset;
        Eigen::Vector3f v3 = p2 + offset;

        vertices.push_back(v0);
        vertices.push_back(v1);
        vertices.push_back(v2);
        vertices.push_back(v3);

        faces.emplace_back(vert_offset, vert_offset + 1, vert_offset + 2);
        faces.emplace_back(vert_offset, vert_offset + 2, vert_offset + 3);
        vert_offset += 4;
    }
}

inline Head::Head(double feature_0, double feature_6)
{
    double width = 1.0 - feature_0 / 2.0;
    double height = 1.5 + feature_0 / 2.0;
    double depth = 0.6;

    auto c = skin_color_from_feature(feature_6);
    auto [v, f] = create_box(width, height, depth);

    set(v, f, c);
}

inline Nose::Nose(double feature_1)
{
    float base = 0.1f;
    float height = 1.0f + feature_1;

    m_color = Eigen::Vector3f(1.0, 0.6, 0.4);
    m_vertices = { { -base, -base, 0 },{ base, -base, 0 },
        { base, base, 0 },{ -base, base, 0 },
        { 0, 0, height } };

    m_faces = { { 0, 1, 4 },{ 1, 2, 4 },{ 2, 3, 4 },{ 3, 0, 4 } };

    translate(0, 0, 0.3);
}

inline Eye::Eye(double feature_5, bool left)
{
    double radius = 0.15;

    Eigen::Vector3f c(0.0, 0.5, 0.8);

    auto [v, f] = create_box(radius, radius, radius);

    set(v, f, c);

    double offset = 0.30 + feature_5 / 5.0;

    offset = left ? -offset : offset;

    translate(offset, 0.3, 0.3);
}

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

//inline QVector<PropWidget*> SwarmOfAvatars3D::create_propwidget_and_get_subwidgets()
//{
//    return PropSwarmOfAvatars3D::create_and_get_subwidgets(this);
//}

