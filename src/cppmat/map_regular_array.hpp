/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_REGULAR_ARRAY_HPP
#define CPPMAT_MAP_REGULAR_ARRAY_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

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

// =================================================================================================
// constructors: map external pointer
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const X *A)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  mData = A;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Map(const X *D)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setMap(D);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U, typename V>
inline
array<X,RANK,I,J,K,L,M,N>::operator std::vector<U> () const
{
  std::vector<U> out(mSize);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = static_cast<U>(mData[i]);

  return out;
}

// =================================================================================================
// modify bounds check
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
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
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return mData[\
    A * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3] +\
    E * mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );
  assert( ( f < nf && f >= -nf ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );
  size_t F = static_cast<size_t>( (nf+(f%nf)) % nf );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3] +\
    E * mStrides[4] +\
    F * mStrides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a) const
{
  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c) const
{
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
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d) const
{
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
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e) const
{
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
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e, T f) const
{
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
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return A * mStrides[0] +\
         B * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return A * mStrides[0] +\
         B * mStrides[1] +\
         C * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );

  return A * mStrides[0] +\
         B * mStrides[1] +\
         C * mStrides[2] +\
         D * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );

  return A * mStrides[0] +\
         B * mStrides[1] +\
         C * mStrides[2] +\
         D * mStrides[3] +\
         E * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );
  assert( ( f < nf && f >= -nf ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );
  size_t F = static_cast<size_t>( (nf+(f%nf)) % nf );

  return A * mStrides[0] +\
         B * mStrides[1] +\
         C * mStrides[2] +\
         D * mStrides[3] +\
         E * mStrides[4] +\
         F * mStrides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a) const
{
  assert( a < mShape[0] );

  return a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return a * mStrides[0] + \
         b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d) const
{
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
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d, T e) const
{
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
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d, T e, T f) const
{
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
// midpoint
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::midpoint() const
{
  // get shape
  std::vector<size_t> mid = shape();

  // check odd-sized
  for ( auto &i : mid )
    if ( i%2 == 0 )
      throw std::domain_error("cppmat::array<X,RANK,I,J,K,L,M,N>::midpoint: Must be odd shaped");

  // midpoint
  for ( auto &i : mid )
    i = (i-1)/2;

  return mid;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::midpoint(size_t axis) const
{
  // get shape
  size_t mid = shape(axis);

  // check odd-sized
  if ( mid%2 == 0 )
    throw std::domain_error("cppmat::array<X,RANK,I,J,K,L,M,N>::midpoint: Must be odd shaped");

  // midpoint
  mid = (mid-1)/2;

  return mid;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X* array<X,RANK,I,J,K,L,M,N>::data() const
{
  return mData;
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::begin() const
{
  return mData;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::end() const
{
  return mData + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

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
auto array<X,RANK,I,J,K,L,M,N>::item(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return begin() +
    A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3] +\
    E * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  assert( ( a < na && a >= -na ) or mPeriodic );
  assert( ( b < nb && b >= -nb ) or mPeriodic );
  assert( ( c < nc && c >= -nc ) or mPeriodic );
  assert( ( d < nd && d >= -nd ) or mPeriodic );
  assert( ( e < ne && e >= -ne ) or mPeriodic );
  assert( ( f < nf && f >= -nf ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );
  size_t D = static_cast<size_t>( (nd+(d%nd)) % nd );
  size_t E = static_cast<size_t>( (ne+(e%ne)) % ne );
  size_t F = static_cast<size_t>( (nf+(f%nf)) % nf );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2] +\
    D * mStrides[3] +\
    E * mStrides[4] +\
    F * mStrides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a) const
{
  assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e, T f) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return begin() +
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
void array<X,RANK,I,J,K,L,M,N>::setMap(const X *D)
{
  mData = D;
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
// bound check
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a, T b) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;
  if ( b >= static_cast<T>(J) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a, T b, T c) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;
  if ( b >= static_cast<T>(J) ) return false;
  if ( c >= static_cast<T>(K) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
    if ( c < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a, T b, T c, T d) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;
  if ( b >= static_cast<T>(J) ) return false;
  if ( c >= static_cast<T>(K) ) return false;
  if ( d >= static_cast<T>(L) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
    if ( c < 0 ) return false;
    if ( d < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a, T b, T c, T d, T e) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;
  if ( b >= static_cast<T>(J) ) return false;
  if ( c >= static_cast<T>(K) ) return false;
  if ( d >= static_cast<T>(L) ) return false;
  if ( e >= static_cast<T>(M) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
    if ( c < 0 ) return false;
    if ( d < 0 ) return false;
    if ( e < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T>
inline
bool array<X,RANK,I,J,K,L,M,N>::inBounds(T a, T b, T c, T d, T e, T f) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(I) ) return false;
  if ( b >= static_cast<T>(J) ) return false;
  if ( c >= static_cast<T>(K) ) return false;
  if ( d >= static_cast<T>(L) ) return false;
  if ( e >= static_cast<T>(M) ) return false;
  if ( f >= static_cast<T>(N) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
    if ( c < 0 ) return false;
    if ( d < 0 ) return false;
    if ( e < 0 ) return false;
    if ( f < 0 ) return false;
  }

  return true;
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
// print operator
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
// equality operators
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
bool operator!= (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
bool operator== (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

