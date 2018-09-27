/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_VECTOR_HPP
#define CPPMAT_FIX_REGULAR_VECTOR_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<typename X, size_t N>
inline
vector<X,N>::vector() : cppmat::tiny::array<X,1,N>()
{
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<typename X, size_t N>
template<typename U, typename V>
inline
vector<X,N>::vector(const cppmat::tiny::array<U,1,N> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<typename X, size_t N>
template<typename U, typename V>
inline
vector<X,N>::vector(const std::vector<U> &A) : cppmat::tiny::array<X,1,N>::Copy(A)
{
}

// =================================================================================================
// constructor: copy from {...}
// =================================================================================================

template<typename X, size_t N>
inline
vector<X,N>::vector(const std::initializer_list<X> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<typename X, size_t N>
inline
vector<X,N>::vector(const cppmat::vector<X> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<typename X, size_t N>
inline
vector<X,N>::vector(const cppmat::view::vector<X,N> &A) : cppmat::tiny::array<X,1,N>(A)
{
}

// =================================================================================================
// finite difference
// =================================================================================================

template<typename X, size_t N>
inline
vector<X,N> vector<X,N>::diff() const
{
  vector<X,N> out;

  std::adjacent_difference(this->begin(), this->end(), out.begin());

  return out;
}

// =================================================================================================
// STL-like behaviour
// =================================================================================================

template<typename X, size_t N>
inline
void vector<X,N>::reserve(size_t n)
{
  Assert( n == N );
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t N>
inline
void vector<X,N>::push_back(const X &value)
{
  Assert( mIstore < N );

  this->mData[mIstore] = value;

  mIstore++;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t N>
inline
void vector<X,N>::clear()
{
  mIstore = 0;

  for ( size_t i = 0 ; i < N ; ++i ) this->mData[i] = static_cast<X>(0);
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

