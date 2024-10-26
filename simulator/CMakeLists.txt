
cmake_minimum_required(VERSION 3.4)
project(simulator_bindings)

find_package(pybind11 REQUIRED)
pybind11_add_module(simulator simulator_bindings.cpp FIELD_3D.cpp PARTICLE_SYSTEM.cpp)
