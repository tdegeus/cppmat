/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_TENSOR2S_HPP
#define CPPMAT_MAP_CARTESIAN_TENSOR2S_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s() : cppmat::view::symmetric::matrix<X,ND,ND>()
{
}

// =================================================================================================
// constructors: map external pointer
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s(const X *A) : cppmat::view::symmetric::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND> tensor2s<X,ND>::Map(const X *D)
{
  tensor2s<X,ND> out;

  out.setMap(D);

  return out;
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
inline
size_t tensor2s<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
