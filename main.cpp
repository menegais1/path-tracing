#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <omp.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "Libs/stb_image.h"
#include "Libs/stb_image_write.h"
#include <glm/gtc/matrix_transform.hpp>

#define byte unsigned char
#define PI 3.145
std::string OUTPUT_PATH = "../Output/";


class Ray {
public:
    glm::dvec3 origin;
    glm::dvec3 direction;

    Ray(glm::dvec3 origin, glm::dvec3 direction) : origin(origin), direction(direction) {}
};

class Material {
public:
    glm::dvec3 color;
    double diffuse;
    glm::dvec3 emission;

    Material(const glm::dvec3 &emission, const glm::dvec3 &color, double diffuse) : color(color), diffuse(diffuse), emission(emission) {}

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

    Sphere(double radius, glm::dvec3 position, Material material) : position(position), radius(radius), material(material) {}
};

class Camera {

public:
    glm::mat4x4 V;

    Camera(glm::vec3 eye, glm::vec3 forward, glm::vec3 up, double fov, double focalLength, int width, int height) : eye(eye), forward(forward), up(up),
                                                                                                                    fov(fov), focalLength(focalLength), width(width), height(height) {
        aspectRatio = width / (float) height;
        auto center = eye + forward;
        V = glm::lookAt(eye, center, up);
    }

    //This function assumes that the camera is positioned at (0,0,0) and pointing in z = -1, so all objects need to be translated from Model to Camera using the View matrix V
    glm::dvec2 rasterToCamera(glm::dvec2 rasterCoord) const {
        // Convert the raster coordinates to [-1,1] range considering the aspect ratio and invert the Y coordinate
        glm::dvec2 screenCoord = glm::vec2(((rasterCoord.x / width) * 2) - 1, (((height - rasterCoord.y) / height) * 2) - 1);
        // Consider the camera field of view and focal length/distance
        glm::dvec2 cameraCoord = screenCoord * tan(focalLength * (fov / 2.0) * (PI / 180));
        return {cameraCoord.x * aspectRatio, cameraCoord.y};
    }

    glm::dvec3 worldToCamera(glm::dvec3 worldCoord) {
        return (V * glm::dvec4(worldCoord, 1));
    }

private:
    double fov;
    double focalLength;
    double aspectRatio;
    glm::dvec3 eye;
    glm::dvec3 up = glm::dvec3(0, 1, 0);
    glm::dvec3 forward;
    int width, height;
};


struct Image {
public:
    int width;
    int height;
    int channels = 3;
    byte *data = new byte[width * height * channels];

    Image(int width, int height) : width(width), height(height) {

    }

    void pixel(int x, int y, glm::dvec3 color) {
        int index = (y * width + x) * channels;
        this->data[index + 0] = (byte) color.x;
        this->data[index + 1] = (byte) color.y;
        this->data[index + 2] = (byte) color.z;

    }

    void write(std::string name) {
        stbi_write_png(name.c_str(), width, height, channels, data, width * channels);
    }
};

bool getClosestHit(Ray r, const std::vector<Hittable *> &sceneObjects, HitInfo &hitInfo) {
    HitInfo curHit{10000};
    bool hit = false;
    for (int i = 0; i < sceneObjects.size(); ++i) {
        Hittable *object = sceneObjects[i];
        if (object->hit(r, curHit)) {
            hit = true;
            if (curHit.t < hitInfo.t) {
                hitInfo = curHit;
            }
        }
    }
    return hit;
}

//glm::dvec3 reflectVector(glm::dvec3 i, glm::dvec3 n) {
//    return i - 2 * glm::dot(i, n) * n;
//}
//
//// https://en.wikipedia.org/wiki/Snell%27s_law
//glm::dvec3 refractVector(glm::dvec3 i, glm::dvec3 n, double n1, double n2) {
//    double nr = n1 / n2;
//    double cosI = -glm::dot(i,n);
//    double totalRefTerm = 1 - (nr * nr) * (1 - (cosI * cosI));
//    if (totalRefTerm < 0)
//        return glm::dvec3(0,0,0);
//    totalRefTerm = sqrt(totalRefTerm);
//    return nr * i + n * (nr * cosI - totalRefTerm);
//}

double uniform_hesmisphere_pdf(glm::dvec2 hemiCoord) {
    return 1.0 / (2 * PI);
//    return sin(hemiCoord[0]) / (2 * PI);
}

glm::dvec3 convert_hemispherical_to_cartesian(glm::dvec2 hemiCoord) {
    double theta = hemiCoord[0];
    double phi = hemiCoord[1];
    return glm::dvec3(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

//returns (theta,phi)
glm::dvec2 uniform_hemisphere_sampling() {
    double u1 = drand48();
    double u2 = drand48();
    return glm::dvec2(acos(u1), u2 * 2 * PI);
}


glm::dvec3 BACKGROUND_COLOR = glm::dvec3(0, 0, 0);

glm::dvec3 traceRay(Ray r, const std::vector<Hittable *> &sceneObjects, int curDepth, int depth = 5) {
    HitInfo info{100000};
    // Diffuse calculation
    if (getClosestHit(r, sceneObjects, info)) {
        if (curDepth > depth) return info.material.emission;
        curDepth++;
        glm::dvec2 hemisphereCoord = uniform_hemisphere_sampling();
        glm::dvec3 disturbance = convert_hemispherical_to_cartesian(hemisphereCoord);
        glm::dvec3 w = glm::normalize(info.normal);
        glm::dvec3 u = glm::normalize(glm::cross(info.normal, abs(info.normal.x) > 0.1 ? glm::dvec3(0, 1, 0) : glm::dvec3(1, 0, 0)));
        glm::dvec3 v = glm::cross(info.normal, u);
        glm::dvec3 sampleDirection = glm::normalize(u * disturbance.x + v * disturbance.y + w * disturbance.z);
        glm::dvec3 diffuse_brdf_coef = info.material.color * (1.0 / PI);
        double cos_term = glm::max(0.0, glm::dot(info.normal, sampleDirection));
        double pdf = uniform_hesmisphere_pdf(hemisphereCoord);
        glm::dvec3 sample_radiance = traceRay(Ray(info.point + sampleDirection * 0.1, sampleDirection), sceneObjects, curDepth);
//        printf("\nSAMPLE RAD: %lf,%lf,%lf",sample_radiance.x,sample_radiance.y,sample_radiance.z);
//        printf("\nPDF: %lf",pdf);
//        printf("\nCOS: %lf",cos_term);
//        printf("\nBRDF: %lf,%lf,%lf",diffuse_brdf_coef.x,diffuse_brdf_coef.y,diffuse_brdf_coef.z);
//        printf("\n---------");
        return info.material.emission + (diffuse_brdf_coef * sample_radiance * cos_term) / pdf;
    }
    return BACKGROUND_COLOR;
}

inline int toInt(double x) {
    return int(pow(glm::clamp(x, 0.0, 1.0), 1 / 2.2) * 255 + .5);
}

int main() {
    std::vector<Hittable *> sceneObjects;
    int width = 1024, height = 768;
    int NUM_SAMPLES = 1000;
    Image output(width, height);

#define SCENE_SMALL_PT

#ifdef SCENE_SMALL_PT
    //smallPt scene description, notice that some things had to me modified to allow running this specific scene, and need to be recoded for running more general scenes, this one is messy.
    Camera camera = Camera(glm::dvec3(50, 52, 295.6), glm::normalize(glm::dvec3(0, -0.042612, -1)), glm::dvec3(0, 1, 0), 29.4213, 1, width, height);
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({1e5 + 1, 40.8, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .25, .25), 1)));//Left
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({-1e5 + 99, 40.8, 81.6}), Material(glm::dvec3(), glm::dvec3(.25, .25, .75), 1)));//Rght
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, 1e5}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75), 1)));//Back
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, -1e5 + 170}), Material(glm::dvec3(), glm::dvec3(0, 0, 0), 1)));//Front Wall.
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 1e5, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75), 1)));//Botm
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, -1e5 + 81.6, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75), 1)));//Top
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({27, 16.5, 47}), Material(glm::dvec3(), glm::dvec3(0.6, 0.34, 1), 1)));//Mirr
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({73, 16.5, 78}), Material(glm::dvec3(), glm::dvec3(1, 0.5, 0.20), 1)));//Glas
    sceneObjects.push_back(new Sphere(600, camera.worldToCamera({50, 681.6 - .27, 81.6}), Material(glm::dvec3(12,12,12), glm::dvec3(), 1))); //Lite
#endif

#ifdef WHITE_FURNACE
    BACKGROUND_COLOR = glm::dvec3(0.5,0.5,0.5);
    Camera camera = Camera(glm::dvec3(0, 0, 0), glm::normalize(glm::dvec3(0, 0, -1)), glm::dvec3(0, 1, 0), 60, 1, width, height);
    sceneObjects.push_back(new Sphere(1, camera.worldToCamera({0, 0, -4}), Material(glm::dvec3(), glm::dvec3(1, 1, 1), 1)));//Left
#endif
    time_t curTime = time(nullptr);

#pragma omp parallel for schedule(dynamic, 1) private(color)
    for (int y = 0; y < output.height; ++y) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", NUM_SAMPLES, 100.0 * y / (height - 1));
        for (int x = 0; x < output.width; ++x) {
            // Create a fixed seed for every line, allowing reproduction from run to run independent of thread order. Taken from SmallPt.
//            unsigned short Xi[3] = {0, 0, static_cast<unsigned short>(y * y * y)};
            glm::dvec3 color = glm::dvec3(0, 0, 0);
            for (int i = 0; i < NUM_SAMPLES; i++) {
                double u = drand48();
                double v = drand48();
                glm::dvec2 cameraCoord = camera.rasterToCamera(glm::dvec2(x + u, y + v));
                glm::dvec3 direction = glm::normalize(glm::dvec3(cameraCoord.x, cameraCoord.y, -1));
                //This shift forward in position is only needed due to smallPt scene description, as one of the spheres overlap the camera at some point.
                Ray r = Ray(glm::dvec3(0, 0, 0) + 140.0 * direction, direction);
                glm::dvec3 radiance = traceRay(r, sceneObjects, 0);
                color = color + radiance;
            }
            color = color / (double) NUM_SAMPLES;

            output.pixel(x, y, glm::vec3(toInt(color.x), toInt(color.y), toInt(color.z)));
        }
    }
    output.write(OUTPUT_PATH + "teste.png");

    fprintf(stderr, "\nRender time: %ld seconds", time(nullptr) - curTime);
    return 0;
}
