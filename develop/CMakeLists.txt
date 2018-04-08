cmake_minimum_required(VERSION 2.8.12)
project(cppmatTest)

# set C++ standard
# - compiler: ... -std=c++14
set(CMAKE_CXX_STANDARD 14)

# option switch on warnings, run : $ cmake .. -DWARNINGS=ON
option(WARNINGS "Show build warnings" ON)
if(WARNINGS)
  if(MSVC)
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
  else()
    # Update if necessary
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  endif()
endif()

# load pkg-config
find_package(PkgConfig)

# load header-only modules using pkg-config
# - cppmat
pkg_check_modules(CPPMAT REQUIRED cppmat)
include_directories(${CPPMAT_INCLUDE_DIRS})
# - eigen3
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

# create executable
add_executable(${PROJECT_NAME}
  matrix.cpp
  matrix2.cpp
  vector.cpp
  tensor.cpp
  tensor2.cpp
  tensor3.cpp
  tiny_matrix2.cpp
  tiny_vector.cpp
)
