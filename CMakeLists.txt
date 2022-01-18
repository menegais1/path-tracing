cmake_minimum_required(VERSION 3.17)
project(PathTracing)

set(CMAKE_CXX_STANDARD 17)

INCLUDE_DIRECTORIES(Libs/glm)
ADD_SUBDIRECTORY(Libs/glm)
add_executable(PathTracing main.cpp Libs/tiny_obj_loader.h Camera.cpp Camera.h Geometry.h Image.cpp Image.h ObjLoader.h)
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
target_link_libraries(PathTracing PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(PathTracing PRIVATE "${OpenMP_CXX_FLAGS}")