/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_ARRAY_CPP
#define CPPMAT_FIX_REGULAR_ARRAY_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// return size without constructing
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::Size()
{
  return I*J*K*L*M*N;
}

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array()
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const array<X,RANK,I,J,K,L,M,N> &A)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  setCopy(A.begin(), A.end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const std::vector<X> &D)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  setCopy(D.begin(), D.end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Random(X lower, X upper)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Arange()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Zero()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Ones()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Constant(X D)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename Iterator>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(Iterator first)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename Iterator>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(Iterator first, Iterator last)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// copy constructor
// =================================================================================================

#ifndef CPPMAT_NOCONVERT
template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline array<X,RANK,I,J,K,L,M,N>::operator cppmat::array<X> () const
{
  return cppmat::array<X>::Copy(shape(), begin(), end());
}
#endif


// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<X> array<X,RANK,I,J,K,L,M,N>::asVector() const
{
  std::vector<X> out(mSize);

  std::copy(begin(), end(), out.begin());

  return out;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( i  <      static_cast<int>(mRank) );
  assert( i >= -1 * static_cast<int>(mRank) );

  // get number of dimensions as integer
  int n = static_cast<int>(mRank);

  // correct periodic index
  i = ( n + (i%n) ) % n;

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  assert( i < mRank );

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::copy(std::begin(mShape), std::begin(mShape)+mRank, shape.begin());

  return shape;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::strides(bool bytes) const
{
  std::vector<size_t> strides(mRank);

  std::copy(std::begin(mStrides), std::begin(mStrides)+mRank, strides.begin());

  if ( bytes )
    for ( size_t i = 0 ; i < mRank ; ++i )
      strides[i] *= sizeof(X);

  return strides;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a)
{
  assert( mRank >= 1 );

  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a) const
{
  assert( mRank >= 1 );

  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b)
{
  assert( mRank >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b) const
{
  assert( mRank >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c)
{
  assert( mRank >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c) const
{
  assert( mRank >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d)
{
  assert( mRank >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mRank >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( mRank >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mRank >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( mRank >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mRank >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5]];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
X& array<X,RANK,I,J,K,L,M,N>::at(Iterator first, Iterator last)
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0     );
  assert( static_cast<size_t>(last-first) <= mRank );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    assert( (*it) < (*shape) );
    // - update index
    idx += (*it) * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return mData[idx];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
const X& array<X,RANK,I,J,K,L,M,N>::at(Iterator first, Iterator last) const
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0     );
  assert( static_cast<size_t>(last-first) <= mRank );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    assert( (*it) < (*shape) );
    // - update index
    idx += (*it) * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return mData[idx];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a) const
{
  assert( mRank >= 1 );

  assert( a < mShape[0] );

  return a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a, size_t b) const
{
  assert( mRank >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return a * mStrides[0] + \
         b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a, size_t b, size_t c) const
{
  assert( mRank >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mRank >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2] + \
         d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mRank >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2] + \
         d * mStrides[3] + \
         e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mRank >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2] + \
         d * mStrides[3] + \
         e * mStrides[4] + \
         f * mStrides[5];
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(mRank);

  // reconstruct
  for ( size_t j = 0 ; j < mRank ; ++j ) {
    idx[j] = (i - i%mStrides[j]) / mStrides[j];
    i -= idx[j] * mStrides[j];
  }

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X* array<X,RANK,I,J,K,L,M,N>::data()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X* array<X,RANK,I,J,K,L,M,N>::data() const
{
  return std::begin(mData);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::begin()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::begin() const
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::end()
{
  return std::begin(mData) + mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::end() const
{
  return std::begin(mData) + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a)
{
  assert( mRank >= 1 );

  assert( a < mShape[0] );

  return begin() + \
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a) const
{
  assert( mRank >= 1 );

  assert( a < mShape[0] );

  return begin() + \
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b)
{
  assert( mRank >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b) const
{
  assert( mRank >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c)
{
  assert( mRank >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c) const
{
  assert( mRank >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d)
{
  assert( mRank >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mRank >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( mRank >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mRank >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( mRank >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mRank >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return begin() + \
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5];
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setRandom(X lower, X upper)
{
  // type of random number distribution
  std::uniform_real_distribution<X> dist(lower, upper);

  // Mersenne Twister: Good quality random number generator
  std::mt19937 rng;

  // Initialize with non-deterministic seeds
  rng.seed(std::random_device{}());

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = dist(rng);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setArange()
{
  std::iota(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// =================================================================================================
// copy to target
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::copyTo(Iterator first, Iterator last) const
{
  assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// =================================================================================================
// sign change
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::operator- () const
{
  array<X,RANK,I,J,K,L,M,N> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::operator+ () const
{
  array<X,RANK,I,J,K,L,M,N> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = mData[i];

  return out;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator*= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator/= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator+= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator-= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::abs() const
{
  array<X,RANK,I,J,K,L,M,N> out;

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = std::abs(mData[i]);

  return out;
}

// =================================================================================================
// norm
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// return the indices that would sort the array
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<size_t,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::argsort(bool ascending) const
{
  return array<size_t,RANK,I,J,K,L,M,N>(cppmat::argsort((*this).asVector(), ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::max() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// =================================================================================================
// mean
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
double array<X,RANK,I,J,K,L,M,N>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(mSize);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
double array<X,RANK,I,J,K,L,M,N>::average(const array<X,RANK,I,J,K,L,M,N> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::where() const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  std::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] ) {
      out[j] = i;
      ++j;
    }
  }

  return out;
}

// =================================================================================================
// find the plain storage indices of all entries equal to some constant
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::where(X D) const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      ++nnz;

  std::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] == D ) {
      out[j] = i;
      ++j;
    }
  }

  return out;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::ostream& operator<<(std::ostream& out, const array<X,RANK,I,J,K,L,M,N>& src)
{
  auto w = out.width();
  auto p = out.precision();

  if ( src.rank() == 1 )
  {
    for ( size_t j = 0 ; j < src.shape(0) ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(j);
      if ( j != src.shape(0)-1 ) out << ", ";
    }

    return out;
  }

  if ( src.rank() == 2 )
  {
    for ( size_t i = 0 ; i < src.shape(0) ; ++i ) {
      for ( size_t j = 0 ; j < src.shape(1) ; ++j ) {
        out << std::setw(w) << std::setprecision(p) << src(i,j);
        if      ( j != src.shape(1)-1 ) out << ", ";
        else if ( i != src.shape(0)-1 ) out << ";" << std::endl;
        else                            out << ";";
      }
    }

    return out;
  }

  out << "cppmat::array[";

  for ( size_t i = 0 ; i < src.rank()-1 ; ++i )
    out << src.shape(i) << ",";

  out << src.shape(src.rank()-1) << "]";

  return out;
}

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  A *= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  A /= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  A += B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  A -= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (array<X,RANK,I,J,K,L,M,N> A, const X &B)
{
  A *= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (array<X,RANK,I,J,K,L,M,N> A, const X &B)
{
  A /= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (array<X,RANK,I,J,K,L,M,N> A, const X &B)
{
  A += B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (array<X,RANK,I,J,K,L,M,N> A, const X &B)
{
  A -= B;

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (const X &A, array<X,RANK,I,J,K,L,M,N> B)
{
  B *= A;

  return B;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (const X &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (const X &A, array<X,RANK,I,J,K,L,M,N> B)
{
  B += A;

  return B;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (const X &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

