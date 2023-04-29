cmake_minimum_required(VERSION 3.16)
project(nvt)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(nvt
        energies.cpp
        forces.cpp
        functions.cpp
        functions.h
        main.cpp
        particle.h
        vector3d.h)
