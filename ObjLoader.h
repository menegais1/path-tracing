//
// Created by Roberto Menegais on 17/01/2022.
//

#ifndef PATHTRACING_OBJLOADER_H
#define PATHTRACING_OBJLOADER_H

#include <iostream>
#include "Libs/tiny_obj_loader.h"
#include "glm/glm.hpp"

/* Position + Normal + UV Vertex */
struct PnuVertexInput {
public:
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 uv;
    glm::vec3 tangent{};

    PnuVertexInput(const glm::vec3 &position, const glm::vec3 &normal, const glm::vec2 &uv) : position(position), normal(normal), uv(uv), tangent(glm::vec3(0, 0, 0)) {};

    PnuVertexInput() : position(glm::vec3(0, 0, 0)), normal(glm::vec3(0, 0, 0)), uv(glm::vec3(0, 0, 0)), tangent(glm::vec3(0, 0, 0)) {};

    bool operator==(const PnuVertexInput &other) const {
        return position == other.position && uv == other.uv && normal == other.normal;
    };
};


namespace std {
    template<>
    struct hash<PnuVertexInput> {
        size_t operator()(PnuVertexInput const &vertex) const {
            return ((hash<glm::dvec3>()(vertex.position) ^
                     (hash<glm::dvec2>()(vertex.uv) << 1)) >> 1) ^
                   (hash<glm::dvec3>()(vertex.normal) << 1);
        }
    };
}







//void LoadObj(const std::string &inputFile, std::vector<PnuVertexInput> &vertices,
//             std::vector<uint32_t> &indices){
//    tinyobj::ObjReader reader;
//    if (!reader.ParseFromFile(inputFile))
//    {
//        if (!reader.Error().empty())
//        {
//            std::cerr << "TinyObjReader: " << reader.Error();
//        }
//        exit(1);
//    }
//
//    if (!reader.Warning().empty())
//    {
//        std::cout << "TinyObjReader: " << reader.Warning();
//    }
//
//    auto &attrib = reader.GetAttrib();
//    auto &shapes = reader.GetShapes();
//    auto &primaryMesh = shapes[0]; /* Fetch the first shape */
//    std::unordered_map<PnuVertexInput, uint32_t> uniqueVertices{};
//    std::cout << "Loading model: " << primaryMesh.name << std::endl;
//
//    for (const auto &index : primaryMesh.mesh.indices)
//    {
//
//        glm::vec3 pos = {attrib.vertices[3 * index.vertex_index + 0],
//                         attrib.vertices[3 * index.vertex_index + 1],
//                         attrib.vertices[3 * index.vertex_index + 2]};
//
//        glm::vec2 texCoord = {attrib.texcoords[2 * index.texcoord_index + 0],
//                              1.0 - attrib.texcoords[2 * index.texcoord_index + 1]};
//
//        glm::vec3 normal = {attrib.normals[3 * index.normal_index + 0],
//                            attrib.normals[3 * index.normal_index + 1],
//                            attrib.normals[3 * index.normal_index + 2]};
//
//        PnuVertexInput vertex = PnuVertexInput(pos, normal, texCoord);
//
//        if (uniqueVertices.count(vertex) == 0) {
//            uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
//            vertices.push_back(vertex);
//        }
//
//        indices.push_back(uniqueVertices[vertex]);
//    }
//    std::cout << primaryMesh.name << " loaded!" << std::endl;
//}

void LoadObj(const std::string &inputFile, std::vector<PnuVertexInput> &vertices,
             std::vector<uint32_t> &indices) {

    std::string inputfile = "cornell_box.obj";
    tinyobj::ObjReaderConfig reader_config;
    reader_config.mtl_search_path = "./"; // Path to material files

    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config)) {
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

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

                // Check if `normal_index` is zero or positive. negative = no normal data
                if (idx.normal_index >= 0) {
                    tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                    tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                    tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
                }

            }
            index_offset += fv;

            // per-face material
            shapes[s].mesh.material_ids[f];
        }
    }
}
#endif //PATHTRACING_OBJLOADER_H
