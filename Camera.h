//
// Created by Roberto Menegais on 17/01/2022.
//

#ifndef PATHTRACING_CAMERA_H
#define PATHTRACING_CAMERA_H


#include "glm/mat4x4.hpp"
#include "glm/ext/matrix_transform.hpp"
#define PI 3.145

class Camera {

public:
    glm::mat4x4 V;

    Camera(glm::vec3 eye, glm::vec3 forward, glm::vec3 up, double fov, double focalLength, int width, int height);

    //This function assumes that the camera is positioned at (0,0,0) and pointing in z = -1, so all objects need to be translated from Model to Camera using the View matrix V
    glm::dvec2 rasterToCamera(glm::dvec2 rasterCoord) const;

    glm::dvec3 worldToCamera(glm::dvec3 worldCoord);
private:
    double fov;
    double focalLength;
    double aspectRatio;
    glm::dvec3 eye;
    glm::dvec3 up = glm::dvec3(0, 1, 0);
    glm::dvec3 forward;
    int width, height;
};


#endif //PATHTRACING_CAMERA_H
