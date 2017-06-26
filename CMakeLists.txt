project(cppmat)

cmake_minimum_required(VERSION 3.0)

include(CMakePackageConfigHelpers)

set(CPPMAT_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
set(CPPMAT_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}" )

set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cppmat" CACHE PATH "Path where cppmatConfig.cmake is installed")

set(PKGCONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/pkgconfig" CACHE PATH "Directory of cppmat.pc")

# Disable pkgconfig only for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(CPPMAT_BUILD_PKGCONFIG "Build pkg-config .pc file for cppmat" ON)
endif()

if(CPPMAT_BUILD_PKGCONFIG)
  configure_file(cppmat.pc.in cppmat.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cppmat.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

configure_package_config_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/cppmatConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cppmatConfig.cmake
  PATH_VARS CPPMAT_INCLUDE_DIR CPPMAT_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
