
# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# define project information
project(cppmat)

# automatically parse the version number
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/cppmat/macros.h" _cppmat_version_header)
string(REGEX MATCH "define[ \t]+CPPMAT_WORLD_VERSION[ \t]+([0-9]+)" _cppmat_world_version_match "${_cppmat_version_header}")
set(CPPMAT_WORLD_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+CPPMAT_MAJOR_VERSION[ \t]+([0-9]+)" _cppmat_major_version_match "${_cppmat_version_header}")
set(CPPMAT_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+CPPMAT_MINOR_VERSION[ \t]+([0-9]+)" _cppmat_minor_version_match "${_cppmat_version_header}")
set(CPPMAT_MINOR_VERSION "${CMAKE_MATCH_1}")
set(CPPMAT_VERSION_NUMBER ${CPPMAT_WORLD_VERSION}.${CPPMAT_MAJOR_VERSION}.${CPPMAT_MINOR_VERSION})

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# set c++ standard -> c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# set paths where 'cppmat' will be installed
# - root
set(CPPMAT_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
# - headers
set(CPPMAT_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
# - headers
set(INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/cppmat"
  CACHE PATH "Path, relative to CMAKE_PREFIX_PATH, where the cppmat header files are installed"
)
# - CMake configuration
set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/cppmat"
  CACHE PATH "Path where cppmatConfig.cmake is installed"
)
# - pkg-config
set(PKGCONFIG_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  CACHE PATH "Path where cppmat.pc is installed"
)

# list all header files
set(CPPMAT_HEADERS
  cppmat/cppmat.h
  cppmat/macros.h
  cppmat/vector.h
  cppmat/matrix.h
  cppmat/matrix2.h
  cppmat/periodic_vector.h
  cppmat/periodic_matrix.h
  cppmat/periodic_matrix2.h
  cppmat/tiny_vector.h
  cppmat/tiny_matrix2.h
  cppmat/tensor.h
  cppmat/tensor2.h
  cppmat/tensor3.h
  cppmat/pybind11_cppmat.h
  cppmat/pybind11_vector.h
  cppmat/pybind11_matrix.h
  cppmat/pybind11_matrix2.h
  cppmat/pybind11_periodic_vector.h
  cppmat/pybind11_periodic_matrix.h
  cppmat/pybind11_periodic_matrix2.h
  cppmat/pybind11_tiny_vector.h
  cppmat/pybind11_tiny_matrix2.h
  cppmat/pybind11_tensor.h
  cppmat/pybind11_tensor2.h
  cppmat/pybind11_tensor3.h
)

# Disable pkg-config for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(CPPMAT_BUILD_PKGCONFIG "Build pkg-config .pc file for cppmat" ON)
endif()

# build pkg-config file -> fill in relevant values
if(CPPMAT_BUILD_PKGCONFIG)
  configure_file(cppmat.pc.in cppmat.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cppmat.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# configure CMake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cppmatConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cppmatConfig.cmake
  PATH_VARS CPPMAT_INCLUDE_DIR CPPMAT_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# - CMake-file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cppmatConfig.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)
# - headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${CPPMAT_HEADERS}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/cppmat
)

# print information to screen
message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in '${CMAKE_INSTALL_PREFIX}'")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake . -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
