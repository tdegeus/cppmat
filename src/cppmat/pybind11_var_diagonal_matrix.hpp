/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_DIAGONAL_MATRIX_PYBIND11_HPP
#define CPPMAT_VAR_DIAGONAL_MATRIX_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::diagonal::matrix <-> NumPy-array
// =================================================================================================

template<class X> struct type_caster<cppmat::diagonal::matrix<X>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::diagonal::matrix<X>, _("cppmat::diagonal::matrix<X>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<X>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<X, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices)
    auto rank = buf.ndim();
    // - check
    if ( rank != 2 ) return false;

    // - shape of the input array
    size_t m = buf.shape()[0];
    size_t n = buf.shape()[1];

    // - all checks passed : create the proper C++ variable
    value = cppmat::diagonal::matrix<X>::CopyDense(m, n, buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::diagonal::matrix<X>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - convert to dense matrix
    cppmat::matrix<X> tmp = src;

    // - create Python variable (all variables are copied)
    py::array a(std::move(tmp.shape()), std::move(tmp.strides(true)), tmp.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================

}} // namespace pybind11::detail

// -------------------------------------------------------------------------------------------------

#endif
