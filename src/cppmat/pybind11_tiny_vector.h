/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_VECTOR_PYBIND11_H
#define CPPMAT_TINY_VECTOR_PYBIND11_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "tiny_vector.h"
#include "macros.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::vector <-> NumPy-array
// =================================================================================================

template <class T, size_t N> struct type_caster<cppmat::tiny::vector<T,N>>
{
public:

  using tinyvec = cppmat::tiny::vector<T,N>;

  PYBIND11_TYPE_CASTER(tinyvec, _("cppmat::tiny::vector<T,N>"));

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

    // - shape : should be exactly N
    if ( buf.shape()[0] != static_cast<ssize_t>(N) ) return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::tiny::vector<T,N>(buf.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::tiny::vector<T,N>& src, py::return_value_policy policy, py::handle parent
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
