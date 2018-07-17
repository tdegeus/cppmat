/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_VECTOR_PYBIND11_HPP
#define CPPMAT_VAR_REGULAR_VECTOR_PYBIND11_HPP

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::vector <-> NumPy-array
// =================================================================================================

template <typename X> struct type_caster<cppmat::vector<X>> :
  list_caster<cppmat::vector<X>, X> { };

// =================================================================================================

}} // namespace pybind11::detail

// -------------------------------------------------------------------------------------------------

#endif
