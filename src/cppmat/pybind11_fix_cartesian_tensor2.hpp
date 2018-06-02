/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR2_PYBIND11_HPP
#define CPPMAT_FIX_CARTESIAN_TENSOR2_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::cartesian::tensor2 <-> NumPy-array
// =================================================================================================

template<class X, size_t ND> struct type_caster<cppmat::tiny::cartesian::tensor2<X,ND>>
{
public:

  using Arr = cppmat::tiny::cartesian::tensor2<X,ND>;

  PYBIND11_TYPE_CASTER(Arr, _("cppmat::tiny::cartesian::tensor2<X,ND>"));

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
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( static_cast<size_t>(buf.shape()[i]) != ND )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::tiny::cartesian::tensor2<X,ND>::Copy(buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::tiny::cartesian::tensor2<X,ND>& src, py::return_value_policy policy, py::handle parent
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
