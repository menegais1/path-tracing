//
// Created by Roberto Menegais on 17/01/2022.
//

#ifndef PATHTRACING_GEOMETRY_H
#define PATHTRACING_GEOMETRY_H


#include "glm/glm.hpp"

class Ray {
public:
    glm::dvec3 origin;
    glm::dvec3 direction;

    Ray(glm::dvec3 origin, glm::dvec3 direction) : origin(origin), direction(direction) {}
};

class Material {
public:
    glm::dvec3 albedo;
    glm::dvec3 emission;

    Material(const glm::dvec3 &emission, const glm::dvec3 &albedo) : albedo(albedo), emission(emission) {}

    Material() = default;
};

class HitInfo {
public:
    double t;
    //This variable indicates if the ray intersected with the inside or outside of an object
    bool insideObject;
    glm::dvec3 normal;
    glm::dvec3 point;
    Material material;
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
//            printf("Dist: %lf", dist);
            info.normal = denom > 0 ? -normal : normal;
            info.material = material;
            return true;
        }

        return false;
    }

    double area() override {
        return 2 * PI * radius;
    }

    double uniform_pdf() {
        return 1.0 / area();
    }

    glm::mat3x3 orthonormal_basis() {
        normal = glm::normalize(normal);
        glm::dvec3 z = glm::normalize(glm::cross(normal, glm::dvec3(-normal.y, normal.x, 0)));
        glm::dvec3 x = glm::cross(z, normal);
        glm::dvec3 y = glm::cross(z, x);

        return {x, y, z};
    }

    void uniform_surface_sampling(glm::dvec3 &point, glm::dvec3 &normal) {
        double u1 = drand48();
        double u2 = drand48();
        double theta = 2 * u2;
        double r = sqrt(u1);
        normal = this->normal;
        point = position + (basis * glm::dvec3(r * cos(theta), 0, r * sin(theta)));
        normal = glm::normalize(normal);
    }

    Disk(glm::dvec3 position, double radius, glm::dvec3 normal, Material material) : position(position), radius(radius),
                                                                                     normal(normal),
                                                                                     material(material) {
        basis = orthonormal_basis();
    }

};

//class Rectangle : public Hittable {
//public:
//    glm::dvec3 p0;
//    glm::dvec3 p1;
//    Material material;
//
//    bool hit(const Ray &r, HitInfo &info) override {
//        // assuming vectors are all normalized
//        float denom = glm::dot(normal, r.direction);
//        if (abs(denom) > 1e-6) {
//            glm::dvec3 p0l0 = p0 - r.origin;
//            info.t = glm::dot(p0l0, normal) / denom;
//            if (info.t < 0) return false;
//            info.point = r.origin + r.direction * info.t;
//            if (info.point.x > max.x || info.point.y > max.y || info.point.z > max.z || info.point.x < min.x ||
//                info.point.y < min.y || info.point.z < min.z)
//                return false;
//            info.normal = denom > 0 ? normal : -normal;
//            info.material = material;
//            return true;
//        }
//
//        return false;
//    }
//
//    void AABB() {
//        min.x = glm::min(p0.x, p1.x);
//        min.y = glm::min(p0.y, p1.y);
//        min.z = glm::min(p0.z, p1.z);
//        max.x = glm::max(p0.x, p1.x);
//        max.y = glm::max(p0.y, p1.y);
//        max.z = glm::max(p0.z, p1.z);
//        auto tmp = glm::dvec3(p1.x, p0.y, p1.z);
//        normal = glm::cross(p1 - p0, tmp - p0);
//    }
//
//    double area() override {
//        throw "error, not implemented";
//    }
//
//    Rectangle(glm::dvec3 p0, glm::dvec3 p1, Material material) : p0(p0), p1(p1), material(material) {
//        AABB();
//    }
//
//private:
//    glm::dvec3 min;
//    glm::dvec3 max;
//    glm::dvec3 normal;
//};
//
//class Triangle : public Hittable {
//public:
//    glm::dvec3 v0, v1, v2;
//
//    Material material;
//
//    bool hit(const Ray &r, HitInfo &info) override {
//        const float EPSILON = 0.0000001;
//        glm::dvec3 vertex0 = v0;
//        glm::dvec3 vertex1 = v1;
//        glm::dvec3 vertex2 = v2;
//        glm::dvec3 edge1, edge2, h, s, q;
//        float a, f, u, v;
//        edge1 = vertex1 - vertex0;
//        edge2 = vertex2 - vertex0;
//        h = glm::cross(r.direction, edge2);
//        a = glm::dot(edge1, h);
//        if (a > -EPSILON && a < EPSILON)
//            return false;    // This ray is parallel to this triangle.
//        f = 1.0 / a;
//        s = r.origin - vertex0;
//        u = f * glm::dot(s, h);
//        if (u < 0.0 || u > 1.0)
//            return false;
//        q = glm::cross(s, edge1);
//        v = f * glm::dot(r.direction, q);
//        if (v < 0.0 || u + v > 1.0)
//            return false;
//        // At this stage we can compute t to find out where the intersection point is on the line.
//        double t = f * glm::dot(edge2, q);
//        if (t > EPSILON) // ray intersection
//        {
//            info.t = t;
//            info.point = r.origin + r.direction * t;
//            info.material = material;
//            info.normal = glm::cross(edge1, edge2);
//            return true;
//        } else // This means that there is a line intersection but not a ray intersection.
//            return false;
//    }
//
//    Triangle(glm::dvec3 v0, glm::dvec3 v1, glm::dvec3 v2, Material material) : v0(v0), v1(v1), v2(v2),
//                                                                               material(material) {}
//
//};


#endif //PATHTRACING_GEOMETRY_H
