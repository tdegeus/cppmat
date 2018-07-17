/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_VECTOR_PYBIND11_HPP
#define CPPMAT_FIX_REGULAR_VECTOR_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::tiny::vector <-> NumPy-array
// =================================================================================================

template <typename X, size_t N> struct type_caster<cppmat::tiny::vector<X,N>> :
  list_caster<cppmat::tiny::vector<X,N>, X> { };

// =================================================================================================

}} // namespace pybind11::detail

// -------------------------------------------------------------------------------------------------

#endif
