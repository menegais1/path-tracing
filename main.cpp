#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <omp.h>
#include <glm/gtx/hash.hpp>

#define TINYOBJLOADER_IMPLEMENTATION

#include "Libs/tiny_obj_loader.h"

#include <unordered_map>
#include "Camera.h"
#include "Geometry.h"
#include "Image.h"
#define byte unsigned char
#define PI 3.145
std::string OUTPUT_PATH = "../Output/";





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
}

double cos_weighted_hemisphere_pdf(glm::dvec2 hemiCoord) {
    return cos(hemiCoord[0]) / PI;
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

glm::dvec2 cos_weighted_hemisphere_sampling() {
    double u1 = drand48();
    double u2 = drand48();
    return glm::dvec2(asin(sqrt(u1)), u2 * 2 * PI);
}


glm::dvec3 BACKGROUND_COLOR = glm::dvec3(0, 0, 0);

glm::dvec3 traceRay(Ray r, const std::vector<Hittable *> &sceneObjects,const std::vector<Hittable *> &lightSources, int curDepth, int depth = 5) {
    HitInfo info{100000};
    // Diffuse calculation
    if (getClosestHit(r, sceneObjects, info)) {
        double rr = drand48();
        glm::dvec3 albedo = info.material.albedo;
        double reflectivity = albedo.x > albedo.y && albedo.x > albedo.z ? albedo.x : albedo.y > albedo.z ? albedo.y : albedo.z;
        if (rr > reflectivity) {
            return info.material.emission;
        }
        glm::dvec2 hemisphereCoord = cos_weighted_hemisphere_sampling();
        glm::dvec3 disturbance = convert_hemispherical_to_cartesian(hemisphereCoord);
        glm::dvec3 w = glm::normalize(info.normal);
        glm::dvec3 u = glm::normalize(glm::cross(info.normal, abs(info.normal.x) > 0.1 ? glm::dvec3(0, 1, 0) : glm::dvec3(1, 0, 0)));
        glm::dvec3 v = glm::cross(info.normal, u);
        glm::dvec3 sampleDirection = glm::normalize(u * disturbance.x + v * disturbance.y + w * disturbance.z);
        glm::dvec3 diffuse_brdf_coef = info.material.albedo * (1.0 / PI);
        double cos_term = glm::max(0.0, glm::dot(info.normal, sampleDirection));
        double pdf = cos_weighted_hemisphere_pdf(hemisphereCoord);
        glm::dvec3 sample_radiance = traceRay(Ray(info.point + sampleDirection * 0.1, sampleDirection), sceneObjects,lightSources, curDepth);
        glm::dvec3 radiance = (info.material.emission + (diffuse_brdf_coef * sample_radiance * cos_term) / pdf);
        return radiance / reflectivity;
    }
    return BACKGROUND_COLOR;
}

inline int toInt(double x) {
    return int(pow(glm::clamp(x, 0.0, 1.0), 1 / 2.2) * 255 + .5);
}

int main() {
    std::vector<Hittable *> sceneObjects;
    std::vector<Hittable *> lightSources;
    int width = 1024, height = 768;
    int NUM_SAMPLES = 100;
    Image output(width, height);

#define SCENE_SMALL_PT

#ifdef SCENE_SMALL_PT
    //smallPt scene description, notice that some things had to me modified to allow running this specific scene, and need to be recoded for running more general scenes, this one is messy.
    Camera camera = Camera(glm::dvec3(50, 52, 295.6), glm::normalize(glm::dvec3(0, -0.042612, -1)), glm::dvec3(0, 1, 0), 29.4213, 1, width, height);
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({1e5 + 1, 40.8, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .25, .25))));//Left
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({-1e5 + 99, 40.8, 81.6}), Material(glm::dvec3(), glm::dvec3(.25, .25, .75))));//Rght
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, 1e5}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Back
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, -1e5 + 170}), Material(glm::dvec3(), glm::dvec3(0, 0, 0))));//Front Wall.
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 1e5, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Botm
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, -1e5 + 81.6, 81.6}), Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Top
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({27, 16.5, 47}), Material(glm::dvec3(), glm::dvec3(0.6, 0.34, 0.99))));//Mirr
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({73, 16.5, 78}), Material(glm::dvec3(), glm::dvec3(0.99, 0.5, 0.20))));//Glas
    sceneObjects.push_back(new Sphere(600, camera.worldToCamera({50, 681.6 - .27, 81.6}), Material(glm::dvec3(12, 12, 12), glm::dvec3(0, 0, 0)))); //Light
#endif

#ifdef PLANE_SCENE
    BACKGROUND_COLOR = glm::dvec3(0.5, 0.5, 0.5);
    //smallPt scene description, notice that some things had to me modified to allow running this specific scene, and need to be recoded for running more general scenes, this one is messy.
    Camera camera = Camera(glm::dvec3(0, 0, 0), glm::normalize(glm::dvec3(-0.1, 0, -1)), glm::dvec3(0, 1, 0), 29.4213, 1, width, height);
//    sceneObjects.push_back(new Sphere(1, camera.worldToCamera({0, 0, -10}), Material(glm::dvec3(), glm::dvec3(.75, .25, .25))));//Left
//    sceneObjects.push_back(new Plane(camera.worldToCamera({-8, 0, 0}), camera.worldToCamera({1, 0, 0.3}), Material(glm::dvec3(), glm::dvec3(.25, .25, .75))));//Left
    sceneObjects.push_back(new Rectangle(camera.worldToCamera({0, -1, -10}), camera.worldToCamera({-2, 1, -10}), Material(glm::dvec3(), glm::dvec3(.25, .25, .75))));//Left
    sceneObjects.push_back(new Triangle(camera.worldToCamera({0, 0, -10}), camera.worldToCamera({2, 0, -10}), camera.worldToCamera({0, 2, -10}), Material(glm::dvec3(), glm::dvec3(.25, .25, .75))));//Left
#endif

#ifdef WHITE_FURNACE
    BACKGROUND_COLOR = glm::dvec3(0.5,0.5,0.5);
    Camera camera = Camera(glm::dvec3(0, 0, 0), glm::normalize(glm::dvec3(0, 0, -1)), glm::dvec3(0, 1, 0), 60, 1, width, height);
    sceneObjects.push_back(new Sphere(1, camera.worldToCamera({0, 0, -4}), Material(glm::dvec3(), glm::dvec3(1, 1, 1))));//Left
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
#ifdef SCENE_SMALL_PT
                Ray r = Ray(glm::dvec3(0, 0, 0) + 140.0 * direction, direction);
#else
                Ray r = Ray(glm::dvec3(0, 0, 0) * direction, direction);
#endif
                glm::dvec3 radiance = traceRay(r, sceneObjects,lightSources, 0);
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
