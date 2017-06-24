cmake_minimum_required(VERSION 2.8.12)
project(tensorlib)

set (CMAKE_CXX_STANDARD 14)

find_package(pybind11 REQUIRED)
pybind11_add_module(tensorlib tensorlib.cpp)
