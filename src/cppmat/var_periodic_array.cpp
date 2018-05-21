/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_PERIODIC_ARRAY_CPP
#define CPPMAT_VAR_PERIODIC_ARRAY_CPP

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
array<X>::array(const std::vector<size_t> &shape) : cppmat::array<X>(shape)
{
  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mShapeI[i] = static_cast<int>(this->mShape[i]);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mStridesI[i] = static_cast<int>(this->mStrides[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>::array(const cppmat::array<X> &A) : cppmat::array<X>(A)
{
  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mShapeI[i] = static_cast<int>(this->mShape[i]);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mStridesI[i] = static_cast<int>(this->mStrides[i]);
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
array<X> array<X>::Random(const std::vector<size_t> &shape, X lower, X upper)
{
  array<X> out(shape);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::Arange(const std::vector<size_t> &shape)
{
  array<X> out(shape);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::Zero(const std::vector<size_t> &shape)
{
  array<X> out(shape);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::Ones(const std::vector<size_t> &shape)
{
  array<X> out(shape);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::Constant(const std::vector<size_t> &shape, X D)
{
  array<X> out(shape);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::Copy(const std::vector<size_t> &shape, const std::vector<X> &D)
{
  array<X> out(shape);

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
array<X> array<X>::Copy(const std::vector<size_t> &shape, Iterator first)
{
  array<X> out(shape);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
array<X> array<X>::Copy(const std::vector<size_t> &shape, Iterator first, Iterator last)
{
  array<X> out(shape);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void array<X>::resize(const std::vector<size_t> &shape)
{
  cppmat::array<X>::resize(shape);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mShapeI[i] = static_cast<int>(this->mShape[i]);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mStridesI[i] = static_cast<int>(this->mStrides[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::reshape(const std::vector<size_t> &shape)
{
  cppmat::array<X>::reshape(shape);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mShapeI[i] = static_cast<int>(this->mShape[i]);

  for ( size_t i = 0 ; i < this->MAX_DIM ; ++i )
    mStridesI[i] = static_cast<int>(this->mStrides[i]);
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline
X& array<X>::operator()(int a)
{
  assert( this->mRank >= 1 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];

  return this->mData[\
    a*mStridesI[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a) const
{
  assert( this->mRank >= 1 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];

  return this->mData[\
    a*mStridesI[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b)
{
  assert( this->mRank >= 2 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a, int b) const
{
  assert( this->mRank >= 2 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b, int c)
{
  assert( this->mRank >= 3 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c) const
{
  assert( this->mRank >= 3 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d)
{
  assert( this->mRank >= 4 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d) const
{
  assert( this->mRank >= 4 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d, int e)
{
  assert( this->mRank >= 5 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d, int e) const
{
  assert( this->mRank >= 5 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d, int e, int f)
{
  assert( this->mRank >= 6 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];
  f = ( mShapeI[5] + (f % mShapeI[5]) ) % mShapeI[5];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4] + \
    f*mStridesI[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d, int e, int f) const
{
  assert( this->mRank >= 6 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];
  f = ( mShapeI[5] + (f % mShapeI[5]) ) % mShapeI[5];

  return this->mData[\
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4] + \
    f*mStridesI[5]];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline
size_t array<X>::compress(int a) const
{
  assert( this->mRank >= 1 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];

  return static_cast<size_t>(
    a*mStridesI[0]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b) const
{
  assert( this->mRank >= 2 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];

  return static_cast<size_t>(
    a*mStridesI[0] + \
    b*mStridesI[1]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b, int c) const
{
  assert( this->mRank >= 3 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];

  return static_cast<size_t>(
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d) const
{
  assert( this->mRank >= 4 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];

  return static_cast<size_t>(
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d, int e) const
{
  assert( this->mRank >= 5 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];

  return static_cast<size_t>(
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d, int e, int f) const
{
  assert( this->mRank >= 6 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];
  f = ( mShapeI[5] + (f % mShapeI[5]) ) % mShapeI[5];

  return static_cast<size_t>(
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4] + \
    f*mStridesI[5]);
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X>
template<class Iterator>
inline
X& array<X>::at(Iterator first, Iterator last)
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0           );
  assert( static_cast<size_t>(last-first) <= this->mRank );

  // iterator to shape and stride
  int *shape   = &this->mShapeI  [0];
  int *stride  = &this->mStridesI[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - current index
    int i = (*it);
    // - correct for periodicity
    i = ( (*shape) + (i % (*shape)) ) % (*shape);
    // - check array index
    assert( i < (*shape) );
    // - update index
    idx += i * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return this->mData[idx];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
const X& array<X>::at(Iterator first, Iterator last) const
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0           );
  assert( static_cast<size_t>(last-first) <= this->mRank );

  // iterator to shape and stride
  int *shape   = &this->mShapeI  [0];
  int *stride  = &this->mStridesI[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - current index
    int i = (*it);
    // - correct for periodicity
    i = ( (*shape) + (i % (*shape)) ) % (*shape);
    // - check array index
    assert( i < (*shape) );
    // - update index
    idx += i * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return this->mData[idx];
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline
auto array<X>::item(int a)
{
  assert( this->mRank >= 1 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];

  return this->begin() + \
    a*mStridesI[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a) const
{
  assert( this->mRank >= 1 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];

  return this->begin() + \
    a*mStridesI[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b)
{
  assert( this->mRank >= 2 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b) const
{
  assert( this->mRank >= 2 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c)
{
  assert( this->mRank >= 3 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c) const
{
  assert( this->mRank >= 3 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d)
{
  assert( this->mRank >= 4 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d) const
{
  assert( this->mRank >= 4 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e)
{
  assert( this->mRank >= 5 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e) const
{
  assert( this->mRank >= 5 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e, int f)
{
  assert( this->mRank >= 6 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];
  f = ( mShapeI[5] + (f % mShapeI[5]) ) % mShapeI[5];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4] + \
    f*mStridesI[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e, int f) const
{
  assert( this->mRank >= 6 );

  a = ( mShapeI[0] + (a % mShapeI[0]) ) % mShapeI[0];
  b = ( mShapeI[1] + (b % mShapeI[1]) ) % mShapeI[1];
  c = ( mShapeI[2] + (c % mShapeI[2]) ) % mShapeI[2];
  d = ( mShapeI[3] + (d % mShapeI[3]) ) % mShapeI[3];
  e = ( mShapeI[4] + (e % mShapeI[4]) ) % mShapeI[4];
  f = ( mShapeI[5] + (f % mShapeI[5]) ) % mShapeI[5];

  return this->begin() + \
    a*mStridesI[0] + \
    b*mStridesI[1] + \
    c*mStridesI[2] + \
    d*mStridesI[3] + \
    e*mStridesI[4] + \
    f*mStridesI[5];
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

