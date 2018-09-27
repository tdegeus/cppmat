/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_ARRAY_PYBIND11_HPP
#define CPPMAT_VAR_REGULAR_ARRAY_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::array <-> NumPy-array
// =================================================================================================

template<class X> struct type_caster<cppmat::array<X>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::array<X>, _("cppmat::array<X>"));

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
    if ( rank < 1 ) return false;

    // - shape of the input array
    std::vector<size_t> shape(rank);
    // - copy
    for ( ssize_t i = 0 ; i < rank ; i++ ) shape[i] = buf.shape()[i];

    // - all checks passed : create the proper C++ variable
    value = cppmat::array<X>::Copy(shape, buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(const cppmat::array<X>& src,
    py::return_value_policy, py::handle)
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================

}} // namespace pybind11::detail

// -------------------------------------------------------------------------------------------------

#endif
