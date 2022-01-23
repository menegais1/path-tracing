//
// Created by Roberto Menegais on 17/01/2022.
//

#ifndef PATHTRACING_GEOMETRY_H
#define PATHTRACING_GEOMETRY_H


#include "glm/glm.hpp"

struct PnuVertexInput {
public:
    glm::dvec3 position;
    glm::dvec3 normal;
    glm::dvec2 uv;
    glm::dvec3 tangent{};

    PnuVertexInput(const glm::dvec3 &position, const glm::dvec3 &normal, const glm::dvec2 &uv) : position(position),
                                                                                                 normal(normal), uv(uv),
                                                                                                 tangent(glm::dvec3(0,
                                                                                                                    0,
                                                                                                                    0)) {};

    PnuVertexInput() : position(glm::dvec3(0, 0, 0)), normal(glm::dvec3(0, 0, 0)), uv(glm::dvec3(0, 0, 0)),
                       tangent(glm::dvec3(0, 0, 0)) {};

    bool operator==(const PnuVertexInput &other) const {
        return position == other.position && uv == other.uv && normal == other.normal;
    };
};

class Ray {
public:
    glm::dvec3 origin;
    glm::dvec3 direction;

    Ray(glm::dvec3 origin, glm::dvec3 direction) : origin(origin), direction(direction) {}

    Ray() = default;
};

enum class MaterialType {
    Default = 1,
    Mirror = 2,
    Glass = 3,
};

class Material {
public:
    glm::dvec3 albedo;
    glm::dvec3 emission;
    MaterialType materialType;
    double n;

    Material(const glm::dvec3 &emission, const glm::dvec3 &albedo, double n = 1.0,
             MaterialType materialType = MaterialType::Default)
            : albedo(albedo),
              emission(emission), n(n),
              materialType(
                      materialType) {}

    Material() = default;
};

class Hittable;

class HitInfo {
public:
    double t;
    //This variable indicates if the ray intersected with the inside or outside of an object
    bool insideObject;
    glm::dvec3 normal;
    glm::dvec3 point;
    Ray incidentRay;
    Material material;
    Hittable *object;
    double refraction = false;
};

class Hittable {
public:
    virtual bool hit(const Ray &r, HitInfo &info) = 0;

    virtual double area() = 0;

    virtual void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) = 0;

    virtual double uniform_pdf() = 0;
};

class Sphere : public Hittable {
public:
    glm::dvec3 position;
    double radius;
    Material material;

    bool hit(const Ray &r, HitInfo &info) override {
        glm::dvec3 vec1 = (r.origin - position);
        double a = glm::dot(r.direction, r.direction);
        double b = glm::dot(r.direction, vec1) * 2.0;
        double c = glm::dot(vec1, vec1) - (radius * radius);

        double delta = (b * b) - (4 * a * c);
        if (delta < 0) return false;
        delta = sqrt(delta);
        double div = (2 * a);
        double r1 = (-b + delta) / div;
        double r2 = (-b - delta) / div;
        if (r1 < 0) r1 = r2;
        if (r2 < 0) r2 = r1;
        if (r2 < 0 && r1 < 0) return false;
        if (r1 <= r2) {
            info.t = r1;
        } else {
            info.t = r2;
        }
        info.point = r.origin + r.direction * info.t;
        info.normal = glm::normalize(info.point - position);
        info.insideObject = false;
        if (glm::dot(info.normal, r.direction) > 0) {
            info.insideObject = true;
            info.normal = -info.normal;
        }
        info.material = material;
        info.object = this;
        info.incidentRay = r;
        return true;
    }

    double area() override {
        return 4 * PI * radius * radius;
    }

    double uniform_pdf() {
        return 1.0 / area();
    }

    void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) {
        double u1 = drand48();
        double u2 = drand48();
        double theta = acos(2 * u1);
        double phi = u2 * 2 * PI;
        normal = glm::dvec3(radius * sin(theta) * cos(phi), radius * sin(theta) * sin(phi), radius * cos(theta));
        point = position + normal;
        normal = glm::normalize(normal);
    }


    Sphere(double radius, glm::dvec3 position, Material material) : position(position), radius(radius),
                                                                    material(material) {}
};

class Disk : public Hittable {
public:
    glm::dvec3 position;
    glm::dvec3 normal;
    double radius;
    Material material;
    glm::dmat3x3 basis;

    bool hit(const Ray &r, HitInfo &info) override {
        // assuming vectors are all normalized
        double denom = glm::dot(normal, r.direction);
        if (glm::abs(denom) > 1e-6) {
            glm::dvec3 p0l0 = position - r.origin;
            info.t = glm::dot(p0l0, normal) / denom;
            if (info.t < 0) return false;
            info.point = r.origin + r.direction * info.t;
            double dist = glm::distance(info.point, position);
            if (dist > radius) return false;
            info.normal = denom > 0 ? -normal : normal;
            info.material = material;
            info.object = this;
            info.incidentRay = r;
            return true;
        }

        return false;
    }

    double area() override {
        return PI * radius * radius;
    }

    double uniform_pdf() {
        return 1.0 / area();
    }

    glm::mat3x3 orthonormal_basis() {
        normal = glm::normalize(normal);
        glm::dvec3 z = glm::normalize(glm::cross(normal, glm::dvec3(-normal.y, normal.x, 0)));
        glm::dvec3 x = glm::normalize(glm::cross(z, normal));
        glm::dvec3 y = glm::normalize(glm::cross(z, x));

        return {x, y, z};
    }

    void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) {
        double u1 = drand48();
        double u2 = drand48();
        double theta = 2 * PI * u2;
        double r = sqrt(u1);
        normal = this->normal;
        auto tmp = (basis * glm::dvec3(r * cos(theta), 0, r * sin(theta)));
        point = position + tmp * radius;
        normal = glm::normalize(normal);
    }

    Disk(glm::dvec3 position, double radius, glm::dvec3 normal, Material material) : position(position), radius(radius),
                                                                                     normal(normal),
                                                                                     material(material) {
        basis = orthonormal_basis();
    }

};

class Triangle : public Hittable {
public:
    PnuVertexInput v0, v1, v2;
    Material material;

    bool hit(const Ray &r, HitInfo &info) override {
        const float EPSILON = 0.0000001;
        glm::dvec3 vertex0 = v0.position;
        glm::dvec3 vertex1 = v1.position;
        glm::dvec3 vertex2 = v2.position;
        glm::dvec3 edge1, edge2, h, s, q;
        float a, f, u, v;
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        h = glm::cross(r.direction, edge2);
        a = glm::dot(edge1, h);
        if (a > -EPSILON && a < EPSILON)
            return false;    // This ray is parallel to this triangle.
        f = 1.0 / a;
        s = r.origin - vertex0;
        u = f * glm::dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = glm::cross(s, edge1);
        v = f * glm::dot(r.direction, q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        double t = f * glm::dot(edge2, q);
        if (t > EPSILON) // ray intersection
        {
            info.t = t;
            info.point = r.origin + r.direction * t;
            info.material = material;
            glm::dvec3 bar = barCoords(info.point);
            info.normal = glm::normalize(bar[0] * v0.normal + bar[1] * v1.normal + bar[2] * v2.normal);
            info.insideObject = false;
            if (glm::dot(info.normal, r.direction) > 0) {
                info.insideObject = true;
                info.normal = -info.normal;
            }
            return true;
        } else // This means that there is a line intersection but not a ray intersection.
            return false;
    }

    glm::dvec3 barCoords(glm::dvec3 p) {
        double a = area();
        double u = glm::length(glm::cross(v0.position - v2.position, p - v2.position)) / 2.0 / a;
        double v = glm::length(glm::cross(v1.position - v0.position, p - v0.position)) / 2.0 / a;
        return {u, v, 1 - u - v};
    }

    double area() override {
        auto a = glm::length(glm::cross(v1.position - v0.position, v2.position - v0.position)) / 2.0;
        return a;
    }

    void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) override {
        double u = drand48();
        double v = drand48();
        double su0 = std::sqrt(u);
        double b0 = 1 - su0;
        double b1 = v * su0;
        point = v0.position * b0 + v1.position * b1 + v2.position * (1.0 - b0 - b1);
        normal = glm::normalize(b0 * v0.normal + b1 * v1.normal + (1.0 - b0 - b1) * v2.normal);
    }

    double uniform_pdf() override {
        return 1 / area();
    }

    Triangle(PnuVertexInput v0, PnuVertexInput v1, PnuVertexInput v2, Material material) : v0(v0), v1(v1), v2(v2),
                                                                                           material(material) {}

};

class Object : public Hittable {
public:
    std::vector<Triangle *> triangles;
    Material material;

    bool hit(const Ray &r, HitInfo &info) override {
        bool hit = false;
        HitInfo curHit{100000};
        for (int i = 0; i < triangles.size(); ++i) {
            if (triangles[i]->hit(r, curHit)) {
                hit = true;
                if (curHit.t < info.t) info = curHit;
            }
        }
        info.material = material;
        info.object = this;
        info.incidentRay = r;
        if (material.materialType == MaterialType::Mirror && info.insideObject) {
            info.insideObject = false;
            info.normal = -info.normal;
        }
        return hit;
    }

    double area() override {
        double sum = 0;
        for (const auto &item: triangles) {
            sum += item->area();
        }
        return sum;
    }

    void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) override {
        int index = (int) (drand48() * triangles.size());
        triangles[index]->uniform_surface_sampling(point, normal);
    }

    double uniform_pdf() override {
        return 1.0 / area();
    }

    Object(std::vector<Triangle *> triangles, Material material) : triangles(triangles), material(material) {

    }
};


#endif //PATHTRACING_GEOMETRY_H
