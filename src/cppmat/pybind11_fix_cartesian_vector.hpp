/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_VECTOR_PYBIND11_HPP
#define CPPMAT_FIX_CARTESIAN_VECTOR_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::cartesian::vector <-> NumPy-array
// =================================================================================================

template <typename X, size_t ND> struct type_caster<cppmat::tiny::cartesian::vector<X,ND>> :
  list_caster<cppmat::tiny::cartesian::vector<X,ND>, X> { };

// =================================================================================================

}} // namespace pybind11::detail

// -------------------------------------------------------------------------------------------------

#endif
