/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_MATRIX_CPP
#define CPPMAT_FIX_REGULAR_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix() : cppmat::tiny::array<X,2,M,N>()
{
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::array<X,2,M,N> &A) : cppmat::tiny::array<X,2,M,N>(A)
{
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::symmetric::matrix<X,M,N> &A) : cppmat::tiny::matrix<X,M,N>()
{
  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      (*this)(i,j) = A(i,j);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::diagonal::matrix<X,M,N> &A) : cppmat::tiny::matrix<X,M,N>()
{
  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      (*this)(i,j) = A(i,j);
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

