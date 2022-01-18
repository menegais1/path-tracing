//
// Created by Roberto Menegais on 17/01/2022.
//

#include "Camera.h"



Camera::Camera(glm::vec3 eye, glm::vec3 forward, glm::vec3 up, double fov, double focalLength, int width, int height) : eye(eye), forward(forward), up(up),
fov(fov), focalLength(focalLength), width(width), height(height) {
    aspectRatio = width / (float) height;
    auto center = eye + forward;
    V = glm::lookAt(eye, center, up);
}

//This function assumes that the camera is positioned at (0,0,0) and pointing in z = -1, so all objects need to be translated from Model to Camera using the View matrix V
glm::dvec2 Camera::rasterToCamera(glm::dvec2 rasterCoord) const {
    // Convert the raster coordinates to [-1,1] range considering the aspect ratio and invert the Y coordinate
    glm::dvec2 screenCoord = glm::vec2(((rasterCoord.x / width) * 2) - 1, (((height - rasterCoord.y) / height) * 2) - 1);
    // Consider the camera field of view and focal length/distance
    glm::dvec2 cameraCoord = screenCoord * tan(focalLength * (fov / 2.0) * (PI / 180));
    return {cameraCoord.x * aspectRatio, cameraCoord.y};
}

glm::dvec3 Camera::worldToCamera(glm::dvec3 worldCoord) {
    return (V * glm::dvec4(worldCoord, 1));
}