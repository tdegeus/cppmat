cmake_minimum_required(VERSION 2.8.0)

project(test)

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -Wsign-compare -Wunused-parameter")
  endif()
endif()

# load pkg-config
find_package(PkgConfig)

# load "eigen3"
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

# create executable
add_executable(${PROJECT_NAME}
  main.cpp
  var_regular_array.cpp
  var_symmetric_matrix.cpp
  var_diagonal_matrix.cpp
  var_misc_matrix.cpp
  var_cartesian_tensor4.cpp
  var_cartesian_tensor2.cpp
  var_cartesian_tensor2s.cpp
  var_cartesian_tensor2d.cpp
  var_cartesian_vector.cpp
  fix_regular_array.cpp
  fix_symmetric_matrix.cpp
  fix_diagonal_matrix.cpp
  fix_misc_matrix.cpp
  fix_cartesian_tensor4.cpp
  fix_cartesian_tensor2.cpp
  fix_cartesian_tensor2s.cpp
  fix_cartesian_tensor2d.cpp
  fix_cartesian_vector.cpp
  fix_cartesian_tensor4_2.cpp
  fix_cartesian_tensor2_2.cpp
  fix_cartesian_tensor2s_2.cpp
  fix_cartesian_tensor2d_2.cpp
  fix_cartesian_vector_2.cpp
  fix_cartesian_tensor4_3.cpp
  fix_cartesian_tensor2_3.cpp
  fix_cartesian_tensor2s_3.cpp
  fix_cartesian_tensor2d_3.cpp
  fix_cartesian_vector_3.cpp
)
