/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_SYMMETRIC_MATRIX_PYBIND11_HPP
#define CPPMAT_FIX_SYMMETRIC_MATRIX_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::symmetric::matrix <-> NumPy-array
// =================================================================================================

template<typename X, size_t M, size_t N> struct type_caster<cppmat::tiny::symmetric::matrix<X,M,N>>
{
public:

  using Arr = cppmat::tiny::symmetric::matrix<X,M,N>;

  PYBIND11_TYPE_CASTER(Arr, _("cppmat::tiny::symmetric::matrix<X,M,N>"));

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

    // - check shape in each direction
    if ( static_cast<size_t>(buf.shape()[0]) != M ) return false;
    if ( static_cast<size_t>(buf.shape()[1]) != N ) return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::tiny::symmetric::matrix<X,M,N>::CopyDense(buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(const cppmat::tiny::symmetric::matrix<X,M,N>& src,
    py::return_value_policy, py::handle)
  {
    // - convert to dense matrix
    cppmat::tiny::matrix<X,M,N> tmp = src;

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
