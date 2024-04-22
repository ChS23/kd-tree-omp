cmake_minimum_required(VERSION 3.28)
project(KDTree)

set(CMAKE_CXX_STANDARD 26)

find_package(OpenMP)

add_executable(KDTree main.cpp
        src/KDNode.h
        src/KDTree.h
        src/KDTree.cpp
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

if (OpenMP_CXX_FOUND)
    target_link_libraries(KDTree PUBLIC OpenMP::OpenMP_CXX)
endif()
