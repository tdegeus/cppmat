/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2D_PYBIND11_HPP
#define CPPMAT_VAR_CARTESIAN_TENSOR2D_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::cartesian::tensor2d <-> NumPy-array
// =================================================================================================

template<typename X> struct type_caster<cppmat::cartesian::tensor2d<X>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::tensor2d<X>, _("cppmat::cartesian::tensor2d<X>"));

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

    // - read number of dimensions (shape in each direction)
    size_t nd = static_cast<size_t>(buf.shape()[0]);
    // - check
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( static_cast<size_t>(buf.shape()[i]) != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::tensor2d<X>::CopyDense(nd, buf.data(), buf.data()+buf.size());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(const cppmat::cartesian::tensor2d<X>& src,
    py::return_value_policy, py::handle)
  {
    // - convert to dense tensor
    cppmat::cartesian::tensor2<X> tmp = src;

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
