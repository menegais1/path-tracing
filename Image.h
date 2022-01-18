//
// Created by caf on 17/01/2022.
//
#pragma
#ifndef PATHTRACING_IMAGE_H
#define PATHTRACING_IMAGE_H

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include <string>
#include "glm/glm.hpp"

#define byte unsigned char

class Image {
public:
    int width;
    int height;
    int channels = 3;
    byte *data = new byte[width * height * channels];

    Image(int width, int height);

    void pixel(int x, int y, glm::dvec3 color) const;

    void write(std::string name);
};

#endif //PATHTRACING_IMAGE_H
