cmake_minimum_required(VERSION 3.17)
project(PathTracing)

set(CMAKE_CXX_STANDARD 17)

ADD_SUBDIRECTORY(Libs/glm)
add_executable(PathTracing main.cpp)
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
target_link_libraries(PathTracing PRIVATE "${OpenMP_CXX_FLAGS}")
target_compile_options(PathTracing PRIVATE "${OpenMP_CXX_FLAGS}")