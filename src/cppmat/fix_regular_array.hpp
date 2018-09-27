/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_ARRAY_HPP
#define CPPMAT_FIX_REGULAR_ARRAY_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// return size without constructing
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::Size()
{
  return I*J*K*L*M*N;
}

// =================================================================================================
// constructors
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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
// constructors: copy from own class (with different type)
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U, typename V>
inline
array<X,RANK,I,J,K,L,M,N>::array(const array<U,RANK,I,J,K,L,M,N> &A)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = static_cast<X>(A[i]);
}

// =================================================================================================
// constructor: copy from {...}
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const std::initializer_list<X> &A)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  this->setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const cppmat::array<X> &A)
{
  Assert( RANK == A.rank() );

  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  Assert( this->shape() == A.shape() );

  setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const cppmat::view::array<X,RANK,I,J,K,L,M,N> &A)
{
  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  setCopy(A.begin(), A.end());
}

// =================================================================================================
// named constructors
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Random(X lower, X upper)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Arange()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Zero()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Ones()
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Constant(X D)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(const std::vector<X> &D)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename Iterator>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(Iterator first)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename Iterator>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(Iterator first, Iterator last)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  Assert( i  <      static_cast<int>(mRank) );
  Assert( i >= -1 * static_cast<int>(mRank) );

  // get number of dimensions as integer
  int n = static_cast<int>(mRank);

  // correct periodic index
  i = ( n + (i%n) ) % n;

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  Assert( i < mRank );

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::copy(std::begin(mShape), std::begin(mShape)+mRank, shape.begin());

  return shape;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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
// get dimensions using a different return type
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
U array<X,RANK,I,J,K,L,M,N>::size() const
{
  return static_cast<U>(size());
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
U array<X,RANK,I,J,K,L,M,N>::rank() const
{
  return static_cast<U>(rank());
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
U array<X,RANK,I,J,K,L,M,N>::shape(int i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
U array<X,RANK,I,J,K,L,M,N>::shape(size_t i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
std::vector<U> array<X,RANK,I,J,K,L,M,N>::shape() const
{
  std::vector<size_t> A = shape();

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename U>
inline
std::vector<U> array<X,RANK,I,J,K,L,M,N>::strides(bool bytes) const
{
  std::vector<size_t> A = strides(bytes);

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator[](size_t i)
{
  Assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator[](size_t i) const
{
  Assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a)
{
  int na = static_cast<int>(mShape[0]);

  Assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return mData[\
    A * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a) const
{
  int na = static_cast<int>(mShape[0]);

  Assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return mData[\
    A * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return mData[\
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e, int f)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );
  Assert( ( f < nf && f >= -nf ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );
  Assert( ( f < nf && f >= -nf ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a)
{
  Assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a) const
{
  Assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e, T f)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );
  Assert( f < mShape[5] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5]];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
const X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e, T f) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );
  Assert( f < mShape[5] );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
X& array<X,RANK,I,J,K,L,M,N>::at(Iterator first, Iterator last)
{
  // check input
  Assert( static_cast<size_t>(last-first)  > 0     );
  Assert( static_cast<size_t>(last-first) <= mRank );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    Assert( (*it) < (*shape) );
    // - update index
    idx += (*it) * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return mData[idx];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
const X& array<X,RANK,I,J,K,L,M,N>::at(Iterator first, Iterator last) const
{
  // check input
  Assert( static_cast<size_t>(last-first)  > 0     );
  Assert( static_cast<size_t>(last-first) <= mRank );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    Assert( (*it) < (*shape) );
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a) const
{
  int na = static_cast<int>(mShape[0]);

  Assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return A * mStrides[0] +\
         B * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return A * mStrides[0] +\
         B * mStrides[1] +\
         C * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );
  Assert( ( f < nf && f >= -nf ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a) const
{
  Assert( a < mShape[0] );

  return a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );

  return a * mStrides[0] + \
         b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2] + \
         d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d, T e) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2] + \
         d * mStrides[3] + \
         e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
size_t array<X,RANK,I,J,K,L,M,N>::compress(T a, T b, T c, T d, T e, T f) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );
  Assert( f < mShape[5] );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
std::vector<size_t> array<X,RANK,I,J,K,L,M,N>::decompress(size_t i) const
{
  // check input
  Assert( i < mSize );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X* array<X,RANK,I,J,K,L,M,N>::data()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
const X* array<X,RANK,I,J,K,L,M,N>::data() const
{
  return std::begin(mData);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::begin()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::begin() const
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::end()
{
  return std::begin(mData) + mSize;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::end() const
{
  return std::begin(mData) + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::index(size_t i)
{
  Assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::index(size_t i) const
{
  Assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a)
{
  int na = static_cast<int>(mShape[0]);

  Assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return begin() +
    A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a) const
{
  int na = static_cast<int>(mShape[0]);

  Assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return begin() +
    A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );
  size_t B = static_cast<size_t>( (nb+(b%nb)) % nb );
  size_t C = static_cast<size_t>( (nc+(c%nc)) % nc );

  return begin() +
    A * mStrides[0] +\
    B * mStrides[1] +\
    C * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e, int f)
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );
  Assert( ( f < nf && f >= -nf ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e, int f) const
{
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  Assert( ( a < na && a >= -na ) or mPeriodic );
  Assert( ( b < nb && b >= -nb ) or mPeriodic );
  Assert( ( c < nc && c >= -nc ) or mPeriodic );
  Assert( ( d < nd && d >= -nd ) or mPeriodic );
  Assert( ( e < ne && e >= -ne ) or mPeriodic );
  Assert( ( f < nf && f >= -nf ) or mPeriodic );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a)
{
  Assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a) const
{
  Assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e, T f)
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );
  Assert( f < mShape[5] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1] + \
    c * mStrides[2] + \
    d * mStrides[3] + \
    e * mStrides[4] + \
    f * mStrides[5];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<typename T, typename S>
inline
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e, T f) const
{
  Assert( a < mShape[0] );
  Assert( b < mShape[1] );
  Assert( c < mShape[2] );
  Assert( d < mShape[3] );
  Assert( e < mShape[4] );
  Assert( f < mShape[5] );

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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setArange()
{
  std::iota(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
void array<X,RANK,I,J,K,L,M,N>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::setCopy(Iterator first, Iterator last)
{
  Assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// =================================================================================================
// copy to target
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
template<class Iterator>
inline
void array<X,RANK,I,J,K,L,M,N>::copyTo(Iterator first, Iterator last) const
{
  Assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// =================================================================================================
// bound check
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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
// sign change
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::operator- () const
{
  array<X,RANK,I,J,K,L,M,N> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator*= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator/= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator+= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator-= (const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator*= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator/= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator+= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>& array<X,RANK,I,J,K,L,M,N>::operator-= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<size_t,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::argsort(bool ascending) const
{
  return array<size_t,RANK,I,J,K,L,M,N>::Copy(cppmat::argsort((*this).asVector(), ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::max() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
X array<X,RANK,I,J,K,L,M,N>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// =================================================================================================
// mean
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
double array<X,RANK,I,J,K,L,M,N>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(mSize);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
double array<X,RANK,I,J,K,L,M,N>::average(const array<X,RANK,I,J,K,L,M,N> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// return array of booleans, based on condition
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::equal(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::not_equal(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::greater(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::greater_equal(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::less(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::less_equal(const X &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::equal(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::not_equal(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::greater(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::greater_equal(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::less(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<int,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::less_equal(const array<X,RANK,I,J,K,L,M,N> &D) const
{
  array<int,RANK,I,J,K,L,M,N> out = array<int,RANK,I,J,K,L,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D[i] )
      out[i] = 1;

  return out;
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::where(size_t index) const
{
  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] ) {
      if ( j == index ) return i;
      ++j;
    }
  }

  throw std::runtime_error("Out-of-bounds");
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
size_t array<X,RANK,I,J,K,L,M,N>::where(int index) const
{
  int nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  Assert( index < nnz && index >= -nnz );

  index = ( nnz + (index%nnz) ) % nnz;

  int j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] ) {
      if ( j == index ) return i;
      ++j;
    }
  }

  throw std::runtime_error("Out-of-bounds");
}

// =================================================================================================
// print operator
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
bool operator!= (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
bool operator== (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator*
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator-
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator*
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const X &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const X &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const X &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator-
(
  const array<X,RANK,I,J,K,L,M,N> &A,
  const X &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator*
(
  const X &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/
(
  const X &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+
(
  const X &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
{
  array<X,RANK,I,J,K,L,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator-
(
  const X &A,
  const array<X,RANK,I,J,K,L,M,N> &B
)
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

