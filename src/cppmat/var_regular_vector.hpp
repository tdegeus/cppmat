/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_VECTOR_HPP
#define CPPMAT_VAR_REGULAR_VECTOR_HPP

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
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<class X>
template<typename U, typename V>
inline
vector<X>::vector(const cppmat::array<U> &A) : cppmat::array<X>(A)
{
  assert( this->mRank == 1 );
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X>
inline
vector<X>::vector(const std::vector<X> &A) : cppmat::array<X>({A.size()})
{
  this->setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<class X>
template<size_t n>
inline
vector<X>::vector(const cppmat::tiny::vector<X,n> &A) : cppmat::array<X>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X>
template<size_t n>
inline
vector<X>::vector(const cppmat::view::vector<X,n> &A) : cppmat::array<X>(A)
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

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(Iterator first, Iterator last)
{
  vector<X> out(last-first);

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
  cppmat::array<X>::resize({n});
}

// =================================================================================================
// extend
// =================================================================================================

template<class X>
inline
void vector<X>::push_back(const X &D)
{
  this->mData.push_back(D);
  this->mShape[0]++;
  this->mSize++;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void vector<X>::append(const cppmat::array<X> &A)
{
  assert( A.rank() == 1 );

  this->mData.insert(this->end(), A.begin(), A.end());

  this->mShape[0] += A.size();
  this->mSize     += A.size();
}

// =================================================================================================
// discrete difference (x1-x0, x2-x1, ...)
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::diff() const
{
  assert( this->mSize > 0 );

  vector<X> out(this->mSize-1);

  for ( size_t i = 0 ; i < this->mSize-1 ; ++i )
    out[i] = this->mData[i+1] - this->mData[i];

  return out;
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
