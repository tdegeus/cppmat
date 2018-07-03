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
// constructors: copy from own class (with different type)
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N>::array(const cppmat::array<X> &A)
{
  assert( RANK == A.rank() );

  mShape[0] = I;  mStrides[0] = J*K*L*M*N;
  mShape[1] = J;  mStrides[1] = K*L*M*N;
  mShape[2] = K;  mStrides[2] = L*M*N;
  mShape[3] = L;  mStrides[3] = M*N;
  mShape[4] = M;  mStrides[4] = N;
  mShape[5] = N;  mStrides[5] = 1;

  assert( this->shape() == A.shape() );

  setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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
inline
array<X,RANK,I,J,K,L,M,N> array<X,RANK,I,J,K,L,M,N>::Copy(const std::vector<X> &D)
{
  array<X,RANK,I,J,K,L,M,N> out;

  out.setCopy(D.begin(), D.end());

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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(int a, int b, int c, int d, int e, int f)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a)
{
  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e)
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
X& array<X,RANK,I,J,K,L,M,N>::operator()(T a, T b, T c, T d, T e, T f)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e)
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
auto array<X,RANK,I,J,K,L,M,N>::item(int a, int b, int c, int d, int e, int f)
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a)
{
  assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b)
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c)
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d)
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e)
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
auto array<X,RANK,I,J,K,L,M,N>::item(T a, T b, T c, T d, T e, T f)
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
  return array<size_t,RANK,I,J,K,L,M,N>::Copy(cppmat::argsort((*this).asVector(), ascending));
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
// return array of booleans, based on condition
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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
// arithmetic operators: external
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
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

