cmake_minimum_required(VERSION 2.8.12)
project(tensorlib)

# set optimization level
# - aggressive optimization
#   compiler: ... -O3
set(CMAKE_BUILD_TYPE Release)
# - switch off assertions (more risky but faster)
#   compiler: ... -DNDEBUG
add_definitions(-DNDEBUG)

# set C++ standard
# - compiler: ... -std=c++14
set(CMAKE_CXX_STANDARD 14)

# find cppmat
# - compiler: ... -I${CPPMAT_INCLUDE_DIRS}
find_package(PkgConfig)
pkg_check_modules(CPPMAT REQUIRED cppmat)
include_directories(${CPPMAT_INCLUDE_DIRS})

# find pybind11
find_package(pybind11 REQUIRED)
pybind11_add_module(tensorlib tensorlib.cpp)
