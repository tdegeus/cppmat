/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX2_PYBIND11_H
#define CPPMAT_TINY_MATRIX2_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::matrix2 <-> NumPy-array
// =================================================================================================

template <class T, size_t M, size_t N> struct type_caster<cppmat::tiny::matrix2<T,M,N>>
{
public:

  using tinymat2 = cppmat::tiny::matrix2<T,M,N>;

  PYBIND11_TYPE_CASTER(tinymat2, _("cppmat::tiny::matrix2<T,M,N>"));

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

    // - rank of the input array (number of indices)
    auto rank = buf.ndim();
    // - check
    if ( rank != 2 ) return false;

    // - shape of the input array: should be exactly (M,N)
    if ( buf.shape()[0] != static_cast<ssize_t>(M) ) return false;
    if ( buf.shape()[1] != static_cast<ssize_t>(N) ) return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::tiny::matrix2<T,M,N>::Copy(buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::tiny::matrix2<T,M,N>& src, py::return_value_policy policy, py::handle parent
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
