/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PYBIND11_H
#define CPPMAT_PYBIND11_H

#include "cppmat.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "pybind11_var_regular_array.hpp"
#include "pybind11_var_regular_matrix.hpp"
#include "pybind11_var_regular_vector.hpp"
#include "pybind11_var_symmetric_matrix.hpp"
#include "pybind11_var_diagonal_matrix.hpp"
#include "pybind11_var_cartesian_tensor4.hpp"
#include "pybind11_var_cartesian_tensor2.hpp"
#include "pybind11_var_cartesian_tensor2s.hpp"
#include "pybind11_var_cartesian_tensor2d.hpp"
#include "pybind11_var_cartesian_vector.hpp"

#include "pybind11_fix_regular_array.hpp"
#include "pybind11_fix_regular_matrix.hpp"
#include "pybind11_fix_regular_vector.hpp"
#include "pybind11_fix_symmetric_matrix.hpp"
#include "pybind11_fix_diagonal_matrix.hpp"
#include "pybind11_fix_cartesian_tensor4.hpp"
#include "pybind11_fix_cartesian_tensor2.hpp"
#include "pybind11_fix_cartesian_tensor2s.hpp"
#include "pybind11_fix_cartesian_tensor2d.hpp"
#include "pybind11_fix_cartesian_vector.hpp"

#endif
