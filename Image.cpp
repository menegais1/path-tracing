//
// Created by caf on 17/01/2022.
//

#include "Image.h"
#include "Libs/stb_image.h"
#include "Libs/stb_image_write.h"
void Image::write(std::string name) {
    stbi_write_png(name.c_str(), width, height, channels, data, width * channels);
}

void Image::pixel(int x, int y, glm::dvec3 color) const {
    int index = (y * width + x) * channels;
    this->data[index + 0] = (byte) color.x;
    this->data[index + 1] = (byte) color.y;
    this->data[index + 2] = (byte) color.z;

}

Image::Image(int width, int height) : width(width), height(height) {

}
