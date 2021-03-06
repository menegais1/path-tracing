#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <omp.h>
#include <glm/gtx/hash.hpp>

#include <unordered_map>
#include "Camera.h"
#include "Geometry.h"
#include "Image.h"
#include "ObjLoader.h"

#define byte unsigned char
#define PI 3.145
std::string OUTPUT_PATH = "../ooutput/";


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

glm::dvec3 reflectVector(glm::dvec3 i, glm::dvec3 n) {
    return glm::normalize(i - 2 * glm::dot(i, n) * n);
}

// https://en.wikipedia.org/wiki/Snell%27s_law
glm::dvec3 refractVector(glm::dvec3 i, glm::dvec3 n, double n1, double n2) {
    double nr = n1 / n2;
    double cosI = -glm::dot(i, n);
    double totalRefTerm = 1 - (nr * nr) * (1 - (cosI * cosI));
    if (totalRefTerm < 0)
        return reflectVector(i, n);
    totalRefTerm = sqrt(totalRefTerm);
    return glm::normalize(nr * i + n * (nr * cosI - totalRefTerm));
}

// https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
double rSchlick2(glm::dvec3 i, glm::dvec3 n, glm::dvec3 t, double n1, double n2) {
    double r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    double cosX = -glm::dot(n, i);
    if (n1 > n2) {
        cosX = -glm::dot(n, t);
    }
    double x = 1 - cosX;
    return r0 + (1 - r0) * x * x * x * x * x;
}


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

double geometricFactor(glm::dvec3 x, glm::dvec3 Nx, glm::dvec3 y, glm::dvec3 Ny) {
    glm::dvec3 dir = glm::normalize(y - x);
    Nx = glm::normalize(Nx);
    Ny = glm::normalize(Ny);
    auto nxCos = glm::max(0.0, glm::dot(Nx, dir));
    auto nyCos = glm::max(0.0, glm::dot(Ny, -dir));
    auto r = glm::distance(x, y);
    return (nxCos * nyCos) / (r * r);
}

glm::dvec3
traceRay(Ray r, const std::vector<Hittable *> &sceneObjects, const std::vector<Hittable *> &lightSources, int depth,
         HitInfo lastHit);

glm::dvec3
direct_radiance(HitInfo info, const std::vector<Hittable *> &sceneObjects,
                const std::vector<Hittable *> &lightSources) {
    double pl = 1.0 / lightSources.size();
    int lightSourceIndex = (int) (drand48() * lightSources.size());
    Hittable *lightSource = lightSources[lightSourceIndex];
    if (info.object == lightSource) return glm::dvec3(0, 0, 0);
    double pyl = lightSource->uniform_pdf();
    glm::dvec3 y, Ny;
    lightSource->uniform_surface_sampling(y, Ny);
    double G = geometricFactor(info.point, info.normal, y, Ny);

    glm::dvec3 direction = glm::normalize(y - info.point);
    HitInfo shadowHit{100000};
    if (getClosestHit(Ray(info.point + direction * 0.1, direction), sceneObjects, shadowHit)) {
        //Light source is visible to point X
        if (glm::distance(shadowHit.point, y) < 0.001) {
            glm::dvec3 diffuse_brdf_coef = info.material.albedo * (1.0 / PI);
            double pdf = pl * pyl;
            glm::dvec3 radiance = (shadowHit.material.emission * diffuse_brdf_coef * G) / pdf;
            return radiance;
        }
        return glm::dvec3(0, 0, 0);
    }

    return BACKGROUND_COLOR;

}

glm::dvec3 indirect_radiance(HitInfo info, const std::vector<Hittable *> &sceneObjects,
                             const std::vector<Hittable *> &lightSources, int depth) {
    double rr = drand48();
    glm::dvec3 albedo = info.material.albedo;
    double reflectivity =
            albedo.x > albedo.y && albedo.x > albedo.z ? albedo.x : albedo.y > albedo.z ? albedo.y : albedo.z;
    if (rr > reflectivity) {
        return glm::dvec3(0, 0, 0);
    }
    glm::dvec3 radiance = glm::dvec3(0);

    if (info.material.materialType == MaterialType::Default) {
        glm::dvec2 hemisphereCoord = cos_weighted_hemisphere_sampling();
        glm::dvec3 disturbance = convert_hemispherical_to_cartesian(hemisphereCoord);
        glm::dvec3 w = glm::normalize(info.normal);
        glm::dvec3 u = glm::normalize(
                glm::cross(info.normal, abs(info.normal.x) > 0.1 ? glm::dvec3(0, 1, 0) : glm::dvec3(1, 0, 0)));
        glm::dvec3 v = glm::cross(info.normal, u);
        glm::dvec3 sampleDirection = glm::normalize(u * disturbance.x + v * disturbance.y + w * disturbance.z);
        glm::dvec3 diffuse_brdf_coef = info.material.albedo * (1.0 / PI);
        double pdf = cos_weighted_hemisphere_pdf(hemisphereCoord);
        double cos_term = glm::max(0.0, glm::dot(info.normal, sampleDirection));
        glm::dvec3 sample_radiance = traceRay(Ray(info.point + sampleDirection * 0.01, sampleDirection), sceneObjects,
                                              lightSources, depth + 1, info);
        radiance = (diffuse_brdf_coef * sample_radiance * cos_term) / pdf;
    } else if (info.material.materialType == MaterialType::Mirror) {
        if (glm::dot(info.incidentRay.direction, info.normal) > 0) return glm::dvec3(0, 0, 0);
        glm::dvec3 sampleDirection = reflectVector(info.incidentRay.direction, info.normal);
        Ray reflectionRay = Ray(info.point + sampleDirection * 0.01, sampleDirection);
        HitInfo lightReflectionHit{10000};
        glm::dvec3 sample_radiance = glm::dvec3(0);
        sample_radiance = traceRay(reflectionRay, sceneObjects,
                                   lightSources, depth + 1, info);

        radiance = sample_radiance;
    } else if (info.material.materialType == MaterialType::Glass) {
        glm::dvec3 sample_radiance = glm::dvec3(0);

        double n1 = 1.0;
        double n2 = info.material.n;
        if (info.insideObject) {
            n1 = info.material.n;
            n2 = 1.0;
        }
        glm::dvec3 refractionDirection = refractVector(info.incidentRay.direction, info.normal, n1, n2);
        double fresnel_reflectance = rSchlick2(info.incidentRay.direction, info.normal, refractionDirection, n1, n2);
        double rr_refraction = drand48();
        if (rr_refraction < fresnel_reflectance) {
            glm::dvec3 sampleDirection = reflectVector(info.incidentRay.direction, info.normal);
            Ray reflectionRay = Ray(info.point + sampleDirection * 0.01, sampleDirection);

            sample_radiance = traceRay(reflectionRay, sceneObjects,
                                       lightSources, depth + 1, info);

            radiance = sample_radiance;
        } else {
            sample_radiance = traceRay(Ray(info.point + refractionDirection * 0.01, refractionDirection), sceneObjects,
                                       lightSources, depth + 1, info);
            radiance = sample_radiance;
        }
    }

    return radiance / reflectivity;

}


glm::dvec3
traceRay(Ray r, const std::vector<Hittable *> &sceneObjects, const std::vector<Hittable *> &lightSources, int depth,
         HitInfo lastHit) {
    HitInfo info{100000};
    // Diffuse calculation
    if (getClosestHit(r, sceneObjects, info)) {
        glm::dvec3 radiance = info.material.emission;
        if (depth > 0 && lastHit.material.materialType == MaterialType::Default) {
            radiance = glm::dvec3(0);
        }

        glm::dvec3 direct = glm::dvec3(0);
        if (info.material.materialType == MaterialType::Default) {
            direct = direct_radiance(info, sceneObjects, lightSources);
        }
        glm::dvec3 indirect = indirect_radiance(info, sceneObjects, lightSources, depth);

        return (radiance + direct + indirect);
    }
    return BACKGROUND_COLOR;
}

//glm::dvec3
//traceRay(Ray r, const std::vector<Hittable *> &sceneObjects, const std::vector<Hittable *> &lightSources, int depth,
//         HitInfo lastHit) {
//    HitInfo info{100000};
//    // Diffuse calculation
//    if (getClosestHit(r, sceneObjects, info)) {
//        double rr = drand48();
//        glm::dvec3 albedo = info.material.albedo;
//        double reflectivity =
//                albedo.x > albedo.y && albedo.x > albedo.z ? albedo.x : albedo.y > albedo.z ? albedo.y : albedo.z;
//        if (rr > reflectivity) {
//            return info.material.emission;
//        }
//        glm::dvec2 hemisphereCoord = cos_weighted_hemisphere_sampling();
//        glm::dvec3 disturbance = convert_hemispherical_to_cartesian(hemisphereCoord);
//        glm::dvec3 w = glm::normalize(info.normal);
//        glm::dvec3 u = glm::normalize(
//                glm::cross(info.normal, abs(info.normal.x) > 0.1 ? glm::dvec3(0, 1, 0) : glm::dvec3(1, 0, 0)));
//        glm::dvec3 v = glm::cross(info.normal, u);
//        glm::dvec3 sampleDirection = glm::normalize(u * disturbance.x + v * disturbance.y + w * disturbance.z);
//        glm::dvec3 diffuse_brdf_coef = info.material.albedo * (1.0 / PI);
//        double cos_term = glm::max(0.0, glm::dot(info.normal, sampleDirection));
//        double pdf = cos_weighted_hemisphere_pdf(hemisphereCoord);
//        glm::dvec3 sample_radiance = traceRay(Ray(info.point + sampleDirection * 0.1, sampleDirection), sceneObjects,
//                                              lightSources,depth,lastHit);
//        glm::dvec3 radiance = (info.material.emission + (diffuse_brdf_coef * sample_radiance * cos_term) / pdf);
//        return radiance / reflectivity;
//    }
//    return BACKGROUND_COLOR;
//}

inline int toInt(double x) {
//    double tone_map = x / (1.0 + x);
    double tone_map = glm::clamp(x, 0.0, 1.0);
//    return int(tone_map * 255 + .5);
    return int(pow(tone_map, 1 / 2.2) * 255 + .5);
}

int main() {
    std::vector<Hittable *> sceneObjects;
    std::vector<Hittable *> lightSources;
    int width = 768, height = 768;
    Image output(width, height);

    std::string scene = "scene_02";
    int NUM_SAMPLES = 8;
    std::string filename = scene + "/" + std::to_string(NUM_SAMPLES) + "spp.png";
#define SCENE_2

#ifdef SCENE_SMALL_PT
    //smallPt scene description, notice that some things had to me modified to allow running this specific scene, and need to be recoded for running more general scenes, this one is messy.
    Camera camera = Camera(glm::dvec3(50, 52, 295.6), glm::normalize(glm::dvec3(0, -0.042612, -1)),
                           glm::dvec3(0, 1, 0),
                           29.4213, 1, width, height);
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({1e5 + 1, 40.8, 81.6}),
                                      Material(glm::dvec3(), glm::dvec3(.75, .25, .25))));//Left
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({-1e5 + 99, 40.8, 81.6}),
                                      Material(glm::dvec3(), glm::dvec3(.25, .25, .75))));//Rght
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, 1e5}),
                                      Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Back
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 40.8, -1e5 + 170}),
                                      Material(glm::dvec3(), glm::dvec3(0, 0, 0))));//Front Wall.
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, 1e5, 81.6}),
                                      Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Botm
    sceneObjects.push_back(new Sphere(1e5, camera.worldToCamera({50, -1e5 + 81.6, 81.6}),
                                      Material(glm::dvec3(), glm::dvec3(.75, .75, .75))));//Top
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({27, 16.5, 47}),
                                      Material(glm::dvec3(),     glm::dvec3(0.6, 0.34, 0.99), 1.0,
                                               MaterialType::Mirror)));//Mirr
    sceneObjects.push_back(new Sphere(16.5, camera.worldToCamera({73, 16.5, 78}),
                                      Material(glm::dvec3(), glm::dvec3(0.99, 0.5, 0.20), 1.5,
                                               MaterialType::Glass)));//Glas
    sceneObjects.push_back(new Sphere(8.5, camera.worldToCamera({53, 16.5, 100}),
                                      Material(glm::dvec3(), glm::dvec3(0.99, 0.5, 0.20))));//Glas


    Disk *lightSource = new Disk(camera.worldToCamera({50, 681.6 - .27 - 599.8, 85}), 17,
                                 camera.normalToCamera({0, -1, 0}),
                                 Material(glm::dvec3(12, 12, 12), glm::dvec3(0, 0, 0)));
    sceneObjects.push_back(lightSource); //Light
    lightSources.push_back(lightSource); //Light
#endif
#define PINK glm::dvec3(0.88, 0.2, 0.6)

#ifdef SCENE_1
    BACKGROUND_COLOR = glm::dvec3(0);
    glm::dvec3 eye = glm::dvec3(5.1425, 4.9361,-14.538);
    glm::dvec3 look_at = glm::dvec3(5.1425, 4.9361,-13.538);
    Camera camera = Camera(eye, look_at - eye, glm::dvec3(0, 1, 0), 40, 1, width, height);
    Object *largeBox = new Object(LoadObj("../scenes/scene_01/cbox_largebox.001.obj", camera),
                                  Material(glm::dvec3(0), glm::dvec3(.5)));
    Object *smallBox = new Object(LoadObj("../scenes/scene_01/cbox_smallbox.001.obj", camera),
                                  Material(glm::dvec3(0), glm::dvec3(.5)));
    Object *backWall = new Object(LoadObj("../scenes/scene_01/cbox_back.001.obj", camera),
                                  Material(glm::dvec3(), glm::dvec3(.4)));
    Object *greenWall = new Object(LoadObj("../scenes/scene_01/cbox_right.001.obj", camera),
                                   Material(glm::dvec3(), glm::dvec3(0,.5,0)));
    Object *redWall = new Object(LoadObj("../scenes/scene_01/cbox_left.obj", camera),
                                 Material(glm::dvec3(), glm::dvec3(.5,0,0)));
    Object *ceiling = new Object(LoadObj("../scenes/scene_01/cbox_ceiling.001.obj", camera),
                                 Material(glm::dvec3(), glm::dvec3(.4)));
    Object *floor = new Object(LoadObj("../scenes/scene_01/cbox_floor.001.obj", camera),
                               Material(glm::dvec3(), glm::dvec3(.4)));
    Object *luminaire1 = new Object(LoadObj("../scenes/scene_01/cbox_luminaire.001.obj", camera, {0,0,0}),
                                    Material(glm::dvec3(.220,.052,.005) *300.0, glm::dvec3(0)));

    sceneObjects.push_back(largeBox);
    sceneObjects.push_back(smallBox);
    sceneObjects.push_back(backWall);
    sceneObjects.push_back(greenWall);
    sceneObjects.push_back(redWall);
    sceneObjects.push_back(ceiling);
    sceneObjects.push_back(floor);
    sceneObjects.push_back(luminaire1);
    lightSources.push_back(luminaire1);
#endif

#ifdef SCENE_2
    BACKGROUND_COLOR = glm::dvec3(0);
    glm::dvec3 eye = glm::dvec3(-32.871, 4.9361, -14.538);
    glm::dvec3 look_at = glm::dvec3(0, 0, -1);
    Camera camera = Camera(eye, {0, 0, 1}, glm::dvec3(0, 1, 0), 40, 1, width, height);
    Object *cylinder = new Object(LoadObj("../scenes/scene_02/Cylinder.obj", camera),
                                  Material(glm::dvec3(0), glm::dvec3(0.202, 0.8, 0.338)));
    Object *icosphere = new Object(LoadObj("../scenes/scene_02/Icosphere.obj", camera),
                                   Material(glm::dvec3(0), glm::dvec3(1), 1.45, MaterialType::Glass));
    Object *plane = new Object(LoadObj("../scenes/scene_02/Plane.obj", camera),
                               Material(glm::dvec3(0), glm::dvec3(1), 1.0, MaterialType::Mirror));
    Material walls = Material(glm::dvec3(), glm::dvec3(0.425, 0.243, 0.800));
//    Material walls = Material(glm::dvec3(), glm::dvec3(0.4));

    Object *backWall = new Object(LoadObj("../scenes/scene_02/cbox_back.obj", camera),
                                  walls);
    Object *greenWall = new Object(LoadObj("../scenes/scene_02/cbox_left.obj", camera),
                                   walls);
    Object *redWall = new Object(LoadObj("../scenes/scene_02/cbox_right.obj", camera),
                                 walls);
    Object *ceiling = new Object(LoadObj("../scenes/scene_02/cbox_ceiling.obj", camera),
                                 walls);
    Object *floor = new Object(LoadObj("../scenes/scene_02/cbox_floor.obj", camera),
                               walls);
    Object *luminaire1 = new Object(LoadObj("../scenes/scene_02/cbox_luminaire.obj", camera, {0, 0, 0}),
                                    Material(glm::dvec3(0.507, 0.241, 0.141) * 300.0, glm::dvec3(0)));

    sceneObjects.push_back(icosphere);
    sceneObjects.push_back(cylinder);
    sceneObjects.push_back(plane);
    sceneObjects.push_back(backWall);
    sceneObjects.push_back(greenWall);
    sceneObjects.push_back(redWall);
    sceneObjects.push_back(ceiling);
    sceneObjects.push_back(floor);
    sceneObjects.push_back(luminaire1);
    lightSources.push_back(luminaire1);
#endif

//Dragon is too heavy for now
#ifdef SCENE_3
    BACKGROUND_COLOR = glm::dvec3(0);
    glm::dvec3 eye = glm::dvec3(-86.566, 12.648, -35.738);
    glm::dvec3 look_at = glm::dvec3(0, 0, -1);
    Camera camera = Camera(eye, {0, 0, 1}, glm::dvec3(0, 1, 0), 40, 1, width, height);
    Object *dragon = new Object(LoadObj("../scenes/scene_03/dragon_vrip_res2.obj", camera),
                                Material(glm::dvec3(0), glm::dvec3(0.801, 0.382, 0.284)));

    Material walls = Material(glm::dvec3(), glm::dvec3(0.247, 0.385, 0.761));

    Object *backWall = new Object(LoadObj("../scenes/scene_03/fundo_infinito.obj", camera),
                                  walls);

    Object *luminaire1 = new Object(LoadObj("../scenes/scene_03/luz_menor.obj", camera, {0, 0, 0}),
                                    Material(glm::dvec3(0.601, 0.943, 1) * 20.0, glm::dvec3(0)));
    Object *luminaire2 = new Object(LoadObj("../scenes/scene_03/luz_maior.obj", camera, {0, 0, 0}),
                                    Material(glm::dvec3(1, 0.862, 0.703) * 30.0, glm::dvec3(0)));

    sceneObjects.push_back(dragon);
    sceneObjects.push_back(backWall);
    sceneObjects.push_back(luminaire1);
    lightSources.push_back(luminaire1);
    sceneObjects.push_back(luminaire2);
    lightSources.push_back(luminaire2);
#endif

#ifdef SCENE_4
    BACKGROUND_COLOR = glm::dvec3(0);
    glm::dvec3 eye = glm::dvec3(-141.28, 4.9361, -14.538);
    glm::dvec3 look_at = glm::dvec3(0, 0, -1);
    Camera camera = Camera(eye, {0, 0, 1}, glm::dvec3(0, 1, 0), 40, 1, width, height);
//    Object *sphere_front = new Object(LoadObj("../scenes/scene_04/esfera_frente.obj", camera),
//                                      Material(glm::dvec3(0), glm::dvec3(1),1.0, MaterialType::Mirror));
//    Object *sphere_back = new Object(LoadObj("../scenes/scene_04/esfera_tras.obj", camera),
//                                     Material(glm::dvec3(0), glm::dvec3(1), 1.45, MaterialType::Glass));
    Sphere *sphere_front = new Sphere(3.79 /2.0, camera.worldToCamera(glm::dvec3(-143.58,1.916,2.1574)),
                                      Material(glm::dvec3(0), glm::dvec3(.9),1.45, MaterialType::Glass));
    Sphere *sphere_back = new Sphere(3.79 /2.0, camera.worldToCamera(glm::dvec3(-139.27,1.916,6.3204)),
                                      Material(glm::dvec3(0), glm::dvec3(.9),1.0, MaterialType::Mirror));
    Object *backWall = new Object(LoadObj("../scenes/scene_04/cbox_back.003.obj", camera),
                                  Material(glm::dvec3(0),glm::dvec3(0.247,0.385,0.761)));
    Object *ceiling = new Object(LoadObj("../scenes/scene_04/cbox_ceiling.003.obj", camera),
                                 Material(glm::dvec3(0),glm::dvec3(1)));
    Object *floor = new Object(LoadObj("../scenes/scene_04/cbox_floor.003.obj", camera),
                               Material(glm::dvec3(0),glm::dvec3(1)));

    Object *greenWall = new Object(LoadObj("../scenes/scene_04/cbox_left.obj", camera),
                                   Material(glm::dvec3(0),glm::dvec3(.9),1.0,MaterialType::Mirror));
    Object *redWall = new Object(LoadObj("../scenes/scene_04/cbox_right.obj", camera),
                                 Material(glm::dvec3(0),glm::dvec3(.9),1.0,MaterialType::Mirror));
    Object *luminaire1 = new Object(LoadObj("../scenes/scene_04/luz_chao.obj", camera, {0, 0, 0}),
                                    Material(glm::dvec3(0.507, 0.241, 0.141) * 50.0, glm::dvec3(0)));
    Object *luminaire2 = new Object(LoadObj("../scenes/scene_04/luz_teto.obj", camera, {0, 0, 0}),
                                    Material(glm::dvec3(0.507, 0.241, 0.141) * 50.0, glm::dvec3(0)));


    sceneObjects.push_back(backWall);
    sceneObjects.push_back(sphere_back);
    sceneObjects.push_back(sphere_front);
    sceneObjects.push_back(greenWall);
    sceneObjects.push_back(redWall);
    sceneObjects.push_back(ceiling);
    sceneObjects.push_back(floor);
    sceneObjects.push_back(luminaire1);
    lightSources.push_back(luminaire1);
    sceneObjects.push_back(luminaire2);
    lightSources.push_back(luminaire2);
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
                glm::dvec3 radiance = traceRay(r, sceneObjects, lightSources, 0, HitInfo{10000});
                color = color + radiance;
            }
            color = color / (double) NUM_SAMPLES;

            output.pixel(x, y, glm::vec3(toInt(color.x), toInt(color.y), toInt(color.z)));
        }
    }
    output.write(OUTPUT_PATH + filename);

    fprintf(stderr, "\nRender time: %ld seconds", time(nullptr) - curTime);
    return 0;
}
