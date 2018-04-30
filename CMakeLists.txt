
# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# project settings
# ----------------

# name
project(cppmat)

# main header file (which the user includes)
set(main_header src/cppmat/cppmat.h)

# all other header files
set(headers
  src/cppmat/stl.h
  src/cppmat/stl.cpp
  src/cppmat/private.h
  src/cppmat/private.cpp
  src/cppmat/regular_array.h
  src/cppmat/regular_array.cpp
  src/cppmat/regular_matrix.h
  src/cppmat/regular_matrix.cpp
  src/cppmat/regular_vector.h
  src/cppmat/regular_vector.cpp
  src/cppmat/periodic_array.h
  src/cppmat/periodic_array.cpp
  src/cppmat/periodic_matrix.h
  src/cppmat/periodic_matrix.cpp
  src/cppmat/periodic_vector.h
  src/cppmat/periodic_vector.cpp
  src/cppmat/tiny_matrix.h
  src/cppmat/tiny_matrix.cpp
  src/cppmat/tiny_vector.h
  src/cppmat/tiny_vector.cpp
  src/cppmat/tensor.h
  src/cppmat/tensor.cpp
  src/cppmat/tensor2.h
  src/cppmat/tensor2.cpp
  src/cppmat/tensor3.h
  src/cppmat/tensor3.cpp
  src/cppmat/view_tiny_matrix.h
  src/cppmat/view_tiny_matrix.cpp
  src/cppmat/view_tiny_vector.h
  src/cppmat/view_tiny_vector.cpp
  src/cppmat/view_tensor2.h
  src/cppmat/view_tensor2.cpp
  src/cppmat/view_tensor3.h
  src/cppmat/view_tensor3.cpp
  src/cppmat/pybind11.h
  src/cppmat/pybind11_regular_array.h
  src/cppmat/pybind11_regular_matrix.h
  src/cppmat/pybind11_regular_vector.h
  src/cppmat/pybind11_periodic_array.h
  src/cppmat/pybind11_periodic_matrix.h
  src/cppmat/pybind11_periodic_vector.h
  src/cppmat/pybind11_tiny_matrix.h
  src/cppmat/pybind11_tiny_vector.h
  src/cppmat/pybind11_tensor.h
  src/cppmat/pybind11_tensor2.h
  src/cppmat/pybind11_tensor3.h
)

# automatically parse the version number
file(READ "${main_header}" version)
string(REGEX MATCH "define[ \t]+CPPMAT_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+CPPMAT_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+CPPMAT_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(CPPMAT_VERSION_NUMBER ${world_version}.${major_version}.${minor_version})

# paths
# -----

set(CPPMAT_ROOT_DIR          "${CMAKE_INSTALL_PREFIX}")
set(CPPMAT_INCLUDE_DIR       "${CMAKE_INSTALL_INCLUDEDIR}")
set(INCLUDE_INSTALL_DIR      "${CMAKE_INSTALL_INCLUDEDIR}")
set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR    "${CMAKE_INSTALL_DATADIR}/pkgconfig")
set(fcmake                   "${PROJECT_NAME}Config.cmake")
set(fpkg                     "${PROJECT_NAME}.pc")

# options
# -------

# configure pkg-config (default: on)
option(PKGCONFIG "Build pkg-config ${fpkg} file" ON)

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${fpkg} file" OFF)
endif()

# C++ standard
# ------------

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# configure CMake
# ---------------

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  PATH_VARS CPPMAT_INCLUDE_DIR CPPMAT_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# -------

# pkg-config
if(PKGCONFIG)
  configure_file(${fpkg}.in ${fpkg} @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fpkg} DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# CMake
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fcmake} DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

# main header file (which the user includes)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${main_header} DESTINATION ${INCLUDE_INSTALL_DIR}/${PROJECT_NAME})

# all other header files
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${headers} DESTINATION ${INCLUDE_INSTALL_DIR}/${PROJECT_NAME})

# print information to screen
# ---------------------------

message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in")
message(STATUS "|   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "|")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake /path/to/${PROJECT_NAME} -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "|")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
