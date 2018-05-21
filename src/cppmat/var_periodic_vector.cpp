/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_PERIODIC_VECTOR_CPP
#define CPPMAT_VAR_PERIODIC_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
vector<X>::vector(size_t n) : cppmat::periodic::array<X>({n})
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const cppmat::periodic::array<X> &A) : cppmat::periodic::array<X>(A)
{
  assert( this->mRank == 1 );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const cppmat::array<X> &A) : cppmat::periodic::array<X>(A)
{
  assert( this->mRank == 1 );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const std::vector<X> &D) : cppmat::periodic::array<X>({D.size()}, D)
{
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::Random(size_t n, X lower, X upper)
{
  vector<X> out(n);

  out.setRandom(lower, upper);

  return out;
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
inline
vector<X> vector<X>::Copy(const std::vector<X> &D)
{
  vector<X> out(D.size());

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Copy(size_t n, const std::vector<X> &D)
{
  vector<X> out(n);

  out.setCopy(D.begin(), D.end());

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
  cppmat::periodic::array<X>::resize({n});
}

// =================================================================================================
// finite difference
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::diff() const
{
  vector<X> out(this->mSize);

  std::adjacent_difference(this->begin(), this->end(), out.begin());

  return out;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

