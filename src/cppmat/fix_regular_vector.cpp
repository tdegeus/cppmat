/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_VECTOR_CPP
#define CPPMAT_FIX_REGULAR_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N>::vector() : cppmat::tiny::array<X,1,N>()
{
}

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N>::vector(const cppmat::tiny::array<X,1,N> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N>::vector(const std::vector<X> &D) : cppmat::tiny::array<X,1,N>::Copy(D)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N>::vector(const cppmat::vector<X> &A)
{
  assert( N == A.size() );

  this->setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N>::vector(const cppmat::view::vector<X,N> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// finite difference
// =================================================================================================

template<class X, size_t N>
inline
vector<X,N> vector<X,N>::diff() const
{
  vector<X,N> out;

  std::adjacent_difference(this->begin(), this->end(), out.begin());

  return out;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

