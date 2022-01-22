//
// Created by Roberto Menegais on 17/01/2022.
//

#ifndef PATHTRACING_OBJLOADER_H
#define PATHTRACING_OBJLOADER_H

#include <iostream>

#define TINYOBJLOADER_IMPLEMENTATION

#include "Libs/tiny_obj_loader.h"
#include "glm/glm.hpp"
#include "Geometry.h"


PnuVertexInput getVertexInfo(const tinyobj::attrib_t &attrib, const std::vector<tinyobj::shape_t> &shapes, size_t s,
                             size_t index_offset) {// access to vertex
    tinyobj::index_t idx = shapes[s].mesh.indices[index_offset];
    glm::dvec3 vertex = glm::dvec3(attrib.vertices[3 * size_t(idx.vertex_index) + 0],
                                   attrib.vertices[3 * size_t(idx.vertex_index) + 1],
                                   attrib.vertices[3 * size_t(idx.vertex_index) + 2]);
    glm::dvec3 normal = glm::dvec3(0);
    // Check if `normal_index` is zero or positive. negative = no normal data
    if (idx.normal_index >= 0) {
        normal = glm::dvec3(attrib.normals[3 * size_t(idx.normal_index) + 0],
                            attrib.normals[3 * size_t(idx.normal_index) + 1],
                            attrib.normals[3 * size_t(idx.normal_index) + 2]);
    }
    return PnuVertexInput(vertex, normal, glm::dvec2(0));
}

std::vector<Triangle *> LoadObj(const std::string &inputFile, Camera camera, glm::dvec3 translate = glm::dvec3(0,0,0)) {

    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./"; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputFile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "TinyObjReader: " << reader.Warning();
    }


    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();
    auto &materials = reader.GetMaterials();


    if (shapes.size() > 1) {
        std::cerr << "This object contains more than one shape, this is not supported";
        exit(-1);
    }
    // Loop over shapes

    std::vector<Triangle *> triangles;
    // Loop over faces(polygon)
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++) {
        size_t fv = size_t(shapes[0].mesh.num_face_vertices[f]);
        if (fv > 3) {
            std::cerr << "Face has more than 3 vertices, this is not supported";
            exit(-1);
        }
        PnuVertexInput v0 = getVertexInfo(attrib, shapes, 0, index_offset);
        PnuVertexInput v1 = getVertexInfo(attrib, shapes, 0, index_offset + 1);
        PnuVertexInput v2 = getVertexInfo(attrib, shapes, 0, index_offset + 2);

        v0.position = camera.worldToCamera(v0.position + translate);
        v1.position = camera.worldToCamera(v1.position + translate);
        v2.position = camera.worldToCamera(v2.position + translate);
        index_offset += fv;

        // per-face material
        shapes[0].mesh.material_ids[f];
        triangles.push_back(new Triangle(v0, v1, v2, Material()));

    }

    glm::dvec3 min = triangles[0]->v0.position;
    glm::dvec3 max = triangles[0]->v0.position;
    for (const auto &item: triangles) {
        min = glm::min(min, item->v0.position);
        min = glm::min(min, item->v1.position);
        min = glm::min(min, item->v2.position);
        max = glm::max(max, item->v2.position);
        max = glm::max(max, item->v2.position);
        max = glm::max(max, item->v2.position);
    }
    auto size = max - min;
    printf("Mesh total size is: %lf %lf %lf", size.x, size.y, size.z);
    return triangles;
}

#endif //PATHTRACING_OBJLOADER_H
