/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_VECTOR_PYBIND11_H
#define CPPMAT_VAR_CARTESIAN_VECTOR_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::cartesian::vector <-> NumPy-array
// =================================================================================================

template<class X> struct type_caster<cppmat::cartesian::vector<X>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::vector<X>, _("cppmat::cartesian::vector<X>"));

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
    if ( rank != 1 ) return false;

    // - read number of dimensions (shape in each direction)
    size_t nd = static_cast<size_t>(buf.shape()[0]);
    // - check
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( static_cast<size_t>(buf.shape()[i]) != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::vector<X>::Copy(nd, buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::vector<X>& src, py::return_value_policy policy, py::handle parent
  )
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