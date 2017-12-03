/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VECTOR_PYBIND11_H
#define CPPMAT_VECTOR_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::periodic::vector <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::periodic::vector<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::periodic::vector<T>, _("cppmat::periodic::vector<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 2
    auto rank = buf.ndim();
    // - check
    if ( rank != 1 ) return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::periodic::vector<T>(buf.shape()[0]);

    // - copy data
    std::copy(buf.data(), buf.data()+value.size(), value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::periodic::vector<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

}} // namespace pybind11::detail

#endif
