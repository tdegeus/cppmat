/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_PERIODIC_MATRIX_CPP
#define CPPMAT_FIX_PERIODIC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t M, size_t N>
CPPMAT_INLINE
matrix<X,M,N>::matrix() : cppmat::tiny::periodic::array<X,2,M,N>()
{
}

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X, size_t M, size_t N>
CPPMAT_INLINE
matrix<X,M,N>::matrix(const cppmat::tiny::array<X,2,M,N> &A) : cppmat::tiny::periodic::array<X,2,M,N>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t M, size_t N>
CPPMAT_INLINE
matrix<X,M,N>::matrix(const cppmat::periodic::matrix<X> &A) : cppmat::tiny::periodic::array<X,2,M,N>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t M, size_t N>
CPPMAT_INLINE
matrix<X,M,N>::matrix(const cppmat::view::periodic::matrix<X,M,N> &A) : cppmat::tiny::periodic::array<X,2,M,N>(A)
{
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

