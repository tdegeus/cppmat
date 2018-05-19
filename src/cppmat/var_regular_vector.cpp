/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_VECTOR_CPP
#define CPPMAT_VAR_REGULAR_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
vector<X>::vector(size_t n) : cppmat::array<X>({n})
{
  N = n;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const cppmat::array<X> &A) : cppmat::array<X>(A)
{
  assert( this->mRank == 1 );

  N = this->mShape[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const std::vector<X> &D) : cppmat::array<X>({D.size()}, D)
{
  N = this->mShape[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Arange(size_t n)
{
  vector<X> out(n);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Zero(size_t n)
{
  vector<X> out(n);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Ones(size_t n)
{
  vector<X> out(n);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Constant(size_t n, X D)
{
  vector<X> out(n);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(size_t n, Iterator first)
{
  vector<X> out(n);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(size_t n, Iterator first, Iterator last)
{
  vector<X> out(n);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void vector<X>::resize(size_t n)
{
  N = n;

  cppmat::array<X>::resize({n});
}

// =================================================================================================
// finite difference
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::diff() const
{
  vector<X> out(N);

  std::adjacent_difference(this->begin(), this->end(), out.begin());

  return out;
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

