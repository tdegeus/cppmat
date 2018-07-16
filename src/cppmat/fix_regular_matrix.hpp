/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_MATRIX_HPP
#define CPPMAT_FIX_REGULAR_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix() : cppmat::tiny::array<X,2,M,N>()
{
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<typename X, size_t M, size_t N>
template<typename U, typename V>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::array<U,2,M,N> &A) : cppmat::tiny::array<X,2,M,N>(A)
{
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::symmetric::matrix<X,M,N> &A) : cppmat::tiny::matrix<X,M,N>()
{
  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      (*this)(i,j) = A(i,j);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::tiny::diagonal::matrix<X,M,N> &A) : cppmat::tiny::matrix<X,M,N>()
{
  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      (*this)(i,j) = A(i,j);
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::matrix<X> &A) : cppmat::tiny::array<X,2,M,N>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::view::matrix<X,M,N> &A) : cppmat::tiny::array<X,2,M,N>(A)
{
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

