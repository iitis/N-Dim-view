#pragma once
#include <Object.h>
#include <Eigen/Dense>

class AvatarPart
{
public:
    std::vector<Eigen::Vector3f> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;
    Eigen::Vector3f m_color;
    Eigen::Vector3f m_position;

    AvatarPart()
    {
        m_vertices = std::vector<Eigen::Vector3f>();
        m_faces = std::vector<Eigen::Vector3i>();
        m_color = { 0,0,0 };
        m_position = { 0, 0, 0 };
    }
    
    void draw(Eigen::Vector3f offset);

    void set(std::vector<Eigen::Vector3f> vertices, std::vector<Eigen::Vector3i> faces, Eigen::Vector3f color)
    {
        m_vertices = vertices;
        m_faces = faces;
        m_color = color;
        m_position = { 0, 0, 0 };
    }

    void translate(float dx, float dy, float dz)
    {
        m_position += Eigen::Vector3f(dx, dy, dz);
    }
};



std::pair<std::vector<Eigen::Vector3f>, std::vector<Eigen::Vector3i>> create_box(float w, float h, float d);
Eigen::Vector3f skin_color_from_feature(float val);


class Head : public AvatarPart
{
public:
    Head(double feature_0, double feature_6)
    {
        double width = 1.0 - feature_0 / 2.0;
        double height = 1.5 + feature_0 / 2.0;
        double depth = 0.6;

        auto c = skin_color_from_feature(feature_6);
        auto [v, f] = create_box(width, height, depth);

        set(v, f, c);
    }
};


class Nose : public AvatarPart
{
public:
    Nose(double feature_1)
    {
        float base = 0.1f;
        float height = 1.0f + feature_1;

        m_color = Eigen::Vector3f(1.0, 0.6, 0.4);
        m_vertices = { {-base, -base, 0}, {base, -base, 0},
            {base, base, 0}, {-base, base, 0},
            {0, 0, height} };

        m_faces = { {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4} };

        translate(0, 0, 0.3);
    }
};

class Eye : public AvatarPart
{
public:
    Eye(double feature_5, bool left = true)
    {
        double radius = 0.15;

        Eigen::Vector3f c(0.0, 0.5, 0.8);

        auto [v, f] = create_box(radius, radius, radius);

        set(v, f, c);

        double offset = 0.30 + feature_5 / 5.0;

        offset = left ? -offset : offset;

        translate(offset, 0.3, 0.3);
    }
};


class Mouth : public AvatarPart
{
public:
    Mouth(double feature_2, double feature_3)
    {
        double width = 0.6 + feature_2 / 2.0;
        double curve = 0.0 + feature_3 / 5.0;
        double thickness = 0.1;
        int segments = 6;

        create_curved_mouth(width, curve, thickness, segments, m_vertices, m_faces);
        m_color = { 1.0, 0.0, 0.0 };

        translate(0, -0.3, 0.31);
    }


    void create_curved_mouth(double width, double curve, double thickness, int segments,
        std::vector<Eigen::Vector3f>& vertices, std::vector<Eigen::Vector3i>& faces)
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
};


class Hair : public AvatarPart
{
public:
    Hair(double feature_4)
    {
        double hair_color = (feature_4 + 1.0) / 2.0;
        Eigen::Vector3f c(hair_color, hair_color, hair_color);
        auto [v, f] = create_box(1.2, 1.0, 0.6);
        set(v, f, c);
        translate(0, 0.4, -0.1);
    }
};


class Avatar3D : public CObject
{
    std::vector<std::shared_ptr<AvatarPart>> m_parts;

    std::vector<Eigen::Vector3f> m_vertices;
    std::vector<Eigen::Vector3i> m_faces;
    std::vector<Eigen::Vector3f> m_colors;

    Eigen::Vector3f m_position;

public:
    Avatar3D(const Eigen::VectorXd& features, const Eigen::Vector3d& pos = Eigen::Vector3d::Zero(4)) : CObject()
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

    void translate(const Eigen::Vector3f& delta) {
        m_position += delta;
    }

    void buildMesh() {
        buildCombinedMesh();
        // VBO creation omitted for brevity
    }

    void draw() const {
        for (const auto& part : m_parts) {
            part->draw(m_position);
        }
    }

    virtual void renderSelf() override
    {
        draw();
    }

    void buildCombinedMesh() {
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

    // VBO-related code would go here
};
