/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_TENSOR2D_CPP
#define CPPMAT_MAP_CARTESIAN_TENSOR2D_CPP

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
tensor2d<X,ND>::tensor2d() : cppmat::view::diagonal::matrix<X,ND,ND>()
{
}

// =================================================================================================
// constructors: map external pointer
// =================================================================================================

template<class X, size_t ND>
inline
tensor2d<X,ND>::tensor2d(const X *A) : cppmat::view::diagonal::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X, size_t ND>
inline
tensor2d<X,ND> tensor2d<X,ND>::Map(const X *D)
{
  tensor2d<X,ND> out;

  out.setMap(D);

  return out;
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
inline
size_t tensor2d<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
