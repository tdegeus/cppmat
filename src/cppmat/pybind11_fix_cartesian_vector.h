/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR2_PYBIND11_H
#define CPPMAT_TENSOR2_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::cartesian::vector <-> NumPy-array
// =================================================================================================

template<class X, size_t ND> struct type_caster<cppmat::tiny::cartesian::vector<X,ND>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::tiny::cartesian::vector<X,ND>, _("cppmat::tiny::cartesian::vector<X,ND>"));

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

    // - check shape in each direction
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( static_cast<size_t>(buf.shape()[i]) != ND )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::tiny::cartesian::vector<X,ND>::Copy(buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::tiny::cartesian::vector<X,ND>& src, py::return_value_policy policy, py::handle parent
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
