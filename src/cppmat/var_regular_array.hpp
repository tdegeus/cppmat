/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_ARRAY_HPP
#define CPPMAT_VAR_REGULAR_ARRAY_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
array<X>::array(const std::vector<size_t> &shape)
{
  resize(shape);
}

// =================================================================================================
// constructors: copy from own class (with different type)
// =================================================================================================

template<class X>
template<typename U, typename V>
inline
array<X>::array(const array<U> &A)
{
  resize(A.shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = static_cast<X>(A[i]);
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<class X>
template<size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X>::array(const cppmat::tiny::array<X,RANK,I,J,K,L,M,N> &A)
{
  resize(A.shape());

  setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X>
template<size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X>::array(const cppmat::view::array<X,RANK,I,J,K,L,M,N> &A)
{
  resize(A.shape());

  setCopy(A.begin(), A.end());
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
// return plain storage as vector
// =================================================================================================

template<class X>
template<typename U, typename V>
inline
array<X>::operator std::vector<U> () const
{
  std::vector<U> out(mSize);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = static_cast<U>(mData[i]);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void array<X>::resize(const std::vector<size_t> &shape)
{
  assert( shape.size() <= MAX_DIM );

  // store old size
  size_t size = mSize;

  // update number of dimensions
  mRank = shape.size();

  // initialize shape/strides in all directions
  std::fill(std::begin(mShape  )+mRank, std::begin(mShape  )+MAX_DIM, 1);
  std::fill(std::begin(mStrides)      , std::begin(mStrides)+MAX_DIM, 1);

  // copy shape from input
  std::copy(shape.begin(), shape.end(), std::begin(mShape));

  // get size
  // - initialize
  mSize = 1;
  // - product of shape in all directions
  for ( size_t i = 0 ; i < mRank ; ++i ) mSize *= shape[i];

  // set storage strides
  for ( size_t i = 0 ; i < MAX_DIM ; ++i )
    for ( size_t j = i+1 ; j < MAX_DIM ; ++j )
      mStrides[i] *= mShape[j];

  // set empty
  if ( shape.size() == 0 ) mSize = 0;

  // allocate data
  if ( mSize != size ) mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::resize(const std::vector<size_t> &shape, const X &D)
{
  assert( shape.size() <= MAX_DIM );

  // store old size
  size_t size = mSize;

  // update number of dimensions
  mRank = shape.size();

  // initialize shape/strides in all directions
  std::fill(std::begin(mShape  )+mRank, std::begin(mShape  )+MAX_DIM, 1);
  std::fill(std::begin(mStrides)      , std::begin(mStrides)+MAX_DIM, 1);

  // copy shape from input
  std::copy(shape.begin(), shape.end(), std::begin(mShape));

  // get size
  // - initialize
  mSize = 1;
  // - product of shape in all directions
  for ( size_t i = 0 ; i < mRank ; ++i ) mSize *= shape[i];

  // set storage strides
  for ( size_t i = 0 ; i < mRank ; ++i )
    for ( size_t j = i+1 ; j < mRank ; ++j )
      mStrides[i] *= mShape[j];

  // set empty
  if ( shape.size() == 0 ) mSize = 0;

  // allocate data
  if ( mSize != size ) mData.resize(mSize, D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::reshape(const std::vector<size_t> &shape)
{
  // check that the size is unchanged
  #ifndef NDEBUG
    size_t n = 1;
    for ( auto &i : shape ) n *= i;
    assert( n == mSize );
  #endif

  // process change
  resize(shape);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::chrank(size_t rank)
{
  // check that all removed dimensions are of shape 1
  #ifndef NDEBUG
    for ( size_t i = rank ; i < MAX_DIM ; ++i ) assert( mShape[i] == 1 );
  #endif

  // update number of dimensions
  mRank = rank;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::ravel()
{
  resize({mSize});
}

// =================================================================================================
// modify bounds check
// =================================================================================================

template<class X>
inline
void array<X>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline
size_t array<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::shape(int i) const
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

template<class X>
inline
size_t array<X>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  assert( i < mRank );

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
std::vector<size_t> array<X>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::copy(std::begin(mShape), std::begin(mShape)+mRank, shape.begin());

  return shape;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
std::vector<size_t> array<X>::strides(bool bytes) const
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

template<class X>
template<typename U>
inline
U array<X>::size() const
{
  return static_cast<U>(size());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U>
inline
U array<X>::rank() const
{
  return static_cast<U>(rank());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U>
inline
U array<X>::shape(int i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U>
inline
U array<X>::shape(size_t i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U>
inline
std::vector<U> array<X>::shape() const
{
  std::vector<size_t> A = shape();

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U>
inline
std::vector<U> array<X>::strides(bool bytes) const
{
  std::vector<size_t> A = strides(bytes);

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X>
inline
X& array<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline
X& array<X>::operator()(int a)
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return mData[\
    A * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X& array<X>::operator()(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return mData[\
    A * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X& array<X>::operator()(int a, int b)
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

template<class X>
inline
const X& array<X>::operator()(int a, int b) const
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

template<class X>
inline
X& array<X>::operator()(int a, int b, int c)
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

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c) const
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

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d)
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

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d) const
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

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d, int e)
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

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d, int e) const
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

template<class X>
inline
X& array<X>::operator()(int a, int b, int c, int d, int e, int f)
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

template<class X>
inline
const X& array<X>::operator()(int a, int b, int c, int d, int e, int f) const
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

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a)
{
  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a) const
{
  assert( a < mShape[0] );

  return mData[\
    a * mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a, T b)
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[\
    a * mStrides[0] + \
    b * mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a, T b, T c)
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

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a, T b, T c) const
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

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a, T b, T c, T d)
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

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a, T b, T c, T d) const
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

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a, T b, T c, T d, T e)
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

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a, T b, T c, T d, T e) const
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

template<class X>
template<typename T, typename S>
inline
X& array<X>::operator()(T a, T b, T c, T d, T e, T f)
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

template<class X>
template<typename T, typename S>
inline
const X& array<X>::operator()(T a, T b, T c, T d, T e, T f) const
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

template<class X>
template<class Iterator>
inline
X& array<X>::at(Iterator first, Iterator last)
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

template<class X>
template<class Iterator>
inline
const X& array<X>::at(Iterator first, Iterator last) const
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

template<class X>
inline
size_t array<X>::compress(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::compress(int a, int b) const
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

template<class X>
inline
size_t array<X>::compress(int a, int b, int c) const
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

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d) const
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

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d, int e) const
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

template<class X>
inline
size_t array<X>::compress(int a, int b, int c, int d, int e, int f) const
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

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a) const
{
  assert( a < mShape[0] );

  return a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return a * mStrides[0] + \
         b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a, T b, T c) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return a * mStrides[0] + \
         b * mStrides[1] + \
         c * mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a, T b, T c, T d) const
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

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a, T b, T c, T d, T e) const
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

template<class X>
template<typename T, typename S>
inline
size_t array<X>::compress(T a, T b, T c, T d, T e, T f) const
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

template<class X>
inline
std::vector<size_t> array<X>::decompress(size_t i) const
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

template<class X>
inline
std::vector<size_t> array<X>::midpoint() const
{
  // get shape
  std::vector<size_t> mid = shape();

  // check odd-sized
  for ( auto &i : mid )
    if ( i%2 == 0 )
      throw std::domain_error("cppmat::array<X>::midpoint: Must be odd shaped");

  // midpoint
  for ( auto &i : mid )
    i = (i-1)/2;

  return mid;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::midpoint(size_t axis) const
{
  // get shape
  size_t mid = shape(axis);

  // check odd-sized
  if ( mid%2 == 0 )
    throw std::domain_error("cppmat::array<X>::midpoint: Must be odd shaped");

  // midpoint
  mid = (mid-1)/2;

  return mid;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline
X* array<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
const X* array<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline
auto array<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline
auto array<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline
auto array<X>::item(int a)
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return begin() +
    A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a) const
{
  int na = static_cast<int>(mShape[0]);

  assert( ( a < na && a >= -na ) or mPeriodic );

  size_t A = static_cast<size_t>( (na+(a%na)) % na );

  return begin() +
    A * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
auto array<X>::item(int a, int b)
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

template<class X>
inline
auto array<X>::item(int a, int b) const
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

template<class X>
inline
auto array<X>::item(int a, int b, int c)
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

template<class X>
inline
auto array<X>::item(int a, int b, int c) const
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d)
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d) const
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e)
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e) const
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e, int f)
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

template<class X>
inline
auto array<X>::item(int a, int b, int c, int d, int e, int f) const
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a)
{
  assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a) const
{
  assert( a < mShape[0] );

  return begin() +
    a * mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b)
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b) const
{
  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() +
    a * mStrides[0] + \
    b * mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c)
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c) const
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d)
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d) const
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d, T e)
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d, T e) const
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d, T e, T f)
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

template<class X>
template<typename T, typename S>
inline
auto array<X>::item(T a, T b, T c, T d, T e, T f) const
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
// slice
// =================================================================================================

template<class X>
inline
array<X> array<X>::slice(
  const std::vector<int> &a, const std::vector<int> &b, const std::vector<int> &c,
  const std::vector<int> &d, const std::vector<int> &e, const std::vector<int> &f
) const
{
  // return empty
  if ( a.size()+b.size()+c.size()+d.size()+e.size()+f.size() == 0 ) return array<X>({0});

  // allocate copies of input lists
  std::vector<size_t> A(a.size());
  std::vector<size_t> B(b.size());
  std::vector<size_t> C(c.size());
  std::vector<size_t> D(d.size());
  std::vector<size_t> E(e.size());
  std::vector<size_t> F(f.size());

  // shape in each direction (as integer)
  int na = static_cast<int>(mShape[0]);
  int nb = static_cast<int>(mShape[1]);
  int nc = static_cast<int>(mShape[2]);
  int nd = static_cast<int>(mShape[3]);
  int ne = static_cast<int>(mShape[4]);
  int nf = static_cast<int>(mShape[5]);

  // check size
  #ifndef NDEBUG
    for ( size_t i = 0 ; i < a.size() ; ++i ) { assert( a[i] < na and a[i] > -na+1 ); }
    for ( size_t i = 0 ; i < b.size() ; ++i ) { assert( b[i] < nb and b[i] > -nb+1 ); }
    for ( size_t i = 0 ; i < c.size() ; ++i ) { assert( c[i] < nc and c[i] > -nc+1 ); }
    for ( size_t i = 0 ; i < d.size() ; ++i ) { assert( d[i] < nd and d[i] > -nd+1 ); }
    for ( size_t i = 0 ; i < e.size() ; ++i ) { assert( e[i] < ne and e[i] > -ne+1 ); }
    for ( size_t i = 0 ; i < f.size() ; ++i ) { assert( f[i] < nf and f[i] > -nf+1 ); }
  #endif

  // copy index from input list: apply periodicity
  for ( size_t i = 0 ; i < a.size() ; ++i ) { A[i] = static_cast<size_t>((na+(a[i]%na))%na); }
  for ( size_t i = 0 ; i < b.size() ; ++i ) { B[i] = static_cast<size_t>((nb+(b[i]%nb))%nb); }
  for ( size_t i = 0 ; i < c.size() ; ++i ) { C[i] = static_cast<size_t>((nc+(c[i]%nc))%nc); }
  for ( size_t i = 0 ; i < d.size() ; ++i ) { D[i] = static_cast<size_t>((nd+(d[i]%nd))%nd); }
  for ( size_t i = 0 ; i < e.size() ; ++i ) { E[i] = static_cast<size_t>((ne+(e[i]%ne))%ne); }
  for ( size_t i = 0 ; i < f.size() ; ++i ) { F[i] = static_cast<size_t>((nf+(f[i]%nf))%nf); }

  // empty list: select all indices
  if ( A.size() == 0 ) { A.resize(mShape[0]); std::iota(A.begin(), A.end(), 0); }
  if ( B.size() == 0 ) { B.resize(mShape[1]); std::iota(B.begin(), B.end(), 0); }
  if ( C.size() == 0 ) { C.resize(mShape[2]); std::iota(C.begin(), C.end(), 0); }
  if ( D.size() == 0 ) { D.resize(mShape[3]); std::iota(D.begin(), D.end(), 0); }
  if ( E.size() == 0 ) { E.resize(mShape[4]); std::iota(E.begin(), E.end(), 0); }
  if ( F.size() == 0 ) { F.resize(mShape[5]); std::iota(F.begin(), F.end(), 0); }

  // shape of the output, without contraction
  // - allocate
  std::vector<size_t> fullshape;
  // - fill
  fullshape.push_back(A.size());
  fullshape.push_back(B.size());
  fullshape.push_back(C.size());
  fullshape.push_back(D.size());
  fullshape.push_back(E.size());
  fullshape.push_back(F.size());

  // allocate output
  array<X> out(fullshape);

  // copy based on selected indices
  for ( size_t i = 0 ; i < A.size() ; ++i )
    for ( size_t j = 0 ; j < B.size() ; ++j )
      for ( size_t k = 0 ; k < C.size() ; ++k )
        for ( size_t l = 0 ; l < D.size() ; ++l )
          for ( size_t m = 0 ; m < E.size() ; ++m )
            for ( size_t n = 0 ; n < F.size() ; ++n )
              out(i,j,k,l,m,n) = (*this)(A[i],B[j],C[k],D[l],E[m],F[n]);

  // shape with contraction
  // - allocate
  std::vector<size_t> shape;
  // - fill
  if ( A.size() > 1 and mRank > 0 ) shape.push_back(A.size());
  if ( B.size() > 1 and mRank > 1 ) shape.push_back(B.size());
  if ( C.size() > 1 and mRank > 2 ) shape.push_back(C.size());
  if ( D.size() > 1 and mRank > 3 ) shape.push_back(D.size());
  if ( E.size() > 1 and mRank > 4 ) shape.push_back(E.size());
  if ( F.size() > 1 and mRank > 5 ) shape.push_back(F.size());
  // - correct for edge case
  if ( shape.size() == 0 ) shape.push_back(1);

  // contract
  out.resize(shape);

  // return output
  return out;
}

// =================================================================================================

template<class X>
template<typename T, typename S>
inline
array<X> array<X>::slice(
  const std::vector<T> &a, const std::vector<T> &b, const std::vector<T> &c,
  const std::vector<T> &d, const std::vector<T> &e, const std::vector<T> &f
) const
{
  // return empty
  if ( a.size()+b.size()+c.size()+d.size()+e.size()+f.size() == 0 ) return array<X>({0});

  // allocate copies of input lists
  std::vector<size_t> A(a.size());
  std::vector<size_t> B(b.size());
  std::vector<size_t> C(c.size());
  std::vector<size_t> D(d.size());
  std::vector<size_t> E(e.size());
  std::vector<size_t> F(f.size());

  // shape in each direction (as integer)
  T na = static_cast<T>(mShape[0]);
  T nb = static_cast<T>(mShape[1]);
  T nc = static_cast<T>(mShape[2]);
  T nd = static_cast<T>(mShape[3]);
  T ne = static_cast<T>(mShape[4]);
  T nf = static_cast<T>(mShape[5]);

  // check size
  #ifndef NDEBUG
    if ( std::numeric_limits<T>::is_signed )
    {
      for ( size_t i = 0 ; i < a.size() ; ++i ) { assert( a[i] < na and a[i] > -na+1 ); }
      for ( size_t i = 0 ; i < b.size() ; ++i ) { assert( b[i] < nb and b[i] > -nb+1 ); }
      for ( size_t i = 0 ; i < c.size() ; ++i ) { assert( c[i] < nc and c[i] > -nc+1 ); }
      for ( size_t i = 0 ; i < d.size() ; ++i ) { assert( d[i] < nd and d[i] > -nd+1 ); }
      for ( size_t i = 0 ; i < e.size() ; ++i ) { assert( e[i] < ne and e[i] > -ne+1 ); }
      for ( size_t i = 0 ; i < f.size() ; ++i ) { assert( f[i] < nf and f[i] > -nf+1 ); }
    }
    else
    {
      for ( size_t i = 0 ; i < a.size() ; ++i ) { assert( a[i] < na and a[i] >= 0 ); }
      for ( size_t i = 0 ; i < b.size() ; ++i ) { assert( b[i] < nb and b[i] >= 0 ); }
      for ( size_t i = 0 ; i < c.size() ; ++i ) { assert( c[i] < nc and c[i] >= 0 ); }
      for ( size_t i = 0 ; i < d.size() ; ++i ) { assert( d[i] < nd and d[i] >= 0 ); }
      for ( size_t i = 0 ; i < e.size() ; ++i ) { assert( e[i] < ne and e[i] >= 0 ); }
      for ( size_t i = 0 ; i < f.size() ; ++i ) { assert( f[i] < nf and f[i] >= 0 ); }
    }
  #endif

  // copy index from input list: apply periodicity
  for ( size_t i = 0 ; i < a.size() ; ++i ) { A[i] = static_cast<size_t>((na+(a[i]%na))%na); }
  for ( size_t i = 0 ; i < b.size() ; ++i ) { B[i] = static_cast<size_t>((nb+(b[i]%nb))%nb); }
  for ( size_t i = 0 ; i < c.size() ; ++i ) { C[i] = static_cast<size_t>((nc+(c[i]%nc))%nc); }
  for ( size_t i = 0 ; i < d.size() ; ++i ) { D[i] = static_cast<size_t>((nd+(d[i]%nd))%nd); }
  for ( size_t i = 0 ; i < e.size() ; ++i ) { E[i] = static_cast<size_t>((ne+(e[i]%ne))%ne); }
  for ( size_t i = 0 ; i < f.size() ; ++i ) { F[i] = static_cast<size_t>((nf+(f[i]%nf))%nf); }

  // empty list: select all indices
  if ( A.size() == 0 ) { A.resize(mShape[0]); std::iota(A.begin(), A.end(), 0); }
  if ( B.size() == 0 ) { B.resize(mShape[1]); std::iota(B.begin(), B.end(), 0); }
  if ( C.size() == 0 ) { C.resize(mShape[2]); std::iota(C.begin(), C.end(), 0); }
  if ( D.size() == 0 ) { D.resize(mShape[3]); std::iota(D.begin(), D.end(), 0); }
  if ( E.size() == 0 ) { E.resize(mShape[4]); std::iota(E.begin(), E.end(), 0); }
  if ( F.size() == 0 ) { F.resize(mShape[5]); std::iota(F.begin(), F.end(), 0); }

  // shape of the output, without contraction
  // - allocate
  std::vector<size_t> fullshape;
  // - fill
  fullshape.push_back(A.size());
  fullshape.push_back(B.size());
  fullshape.push_back(C.size());
  fullshape.push_back(D.size());
  fullshape.push_back(E.size());
  fullshape.push_back(F.size());

  // allocate output
  array<X> out(fullshape);

  // copy based on selected indices
  for ( size_t i = 0 ; i < A.size() ; ++i )
    for ( size_t j = 0 ; j < B.size() ; ++j )
      for ( size_t k = 0 ; k < C.size() ; ++k )
        for ( size_t l = 0 ; l < D.size() ; ++l )
          for ( size_t m = 0 ; m < E.size() ; ++m )
            for ( size_t n = 0 ; n < F.size() ; ++n )
              out(i,j,k,l,m,n) = (*this)(A[i],B[j],C[k],D[l],E[m],F[n]);

  // shape with contraction
  // - allocate
  std::vector<size_t> shape;
  // - fill
  if ( A.size() > 1 and mRank > 0 ) shape.push_back(A.size());
  if ( B.size() > 1 and mRank > 1 ) shape.push_back(B.size());
  if ( C.size() > 1 and mRank > 2 ) shape.push_back(C.size());
  if ( D.size() > 1 and mRank > 3 ) shape.push_back(D.size());
  if ( E.size() > 1 and mRank > 4 ) shape.push_back(E.size());
  if ( F.size() > 1 and mRank > 5 ) shape.push_back(F.size());
  // - correct for edge case
  if ( shape.size() == 0 ) shape.push_back(1);

  // contract
  out.resize(shape);

  // return output
  return out;
}

// =================================================================================================
// return padded array
// =================================================================================================

template<class X>
inline
array<X> array<X>::pad(const std::vector<size_t> &pad_width, X D)
{
  assert( pad_width.size() == mRank );

  // output array
  // - allocate shape
  std::vector<size_t> shape(mRank);
  // - fill shape
  for ( size_t i = 0 ; i < mRank ; ++i )
    shape[i] = mShape[i] + 2*pad_width[i];
  // - allocate array
  array<X> out = array<X>::Constant(shape, D);

  // pad size
  // - allocate
  std::vector<size_t> pad(MAX_DIM, 0);
  // - fill
  for ( size_t i = 0 ; i < mRank ; ++i )
    pad[i] = pad_width[i];

  // place current array in output
  for ( size_t i = 0 ; i < mShape[0] ; ++i )
    for ( size_t j = 0 ; j < mShape[1] ; ++j )
      for ( size_t k = 0 ; k < mShape[2] ; ++k )
        for ( size_t l = 0 ; l < mShape[3] ; ++l )
          for ( size_t m = 0 ; m < mShape[4] ; ++m )
            for ( size_t n = 0 ; n < mShape[5] ; ++n )
              out(i+pad[0],j+pad[1],k+pad[2],l+pad[3],m+pad[4],n+pad[5]) = (*this)(i,j,k,l,m,n);

  return out;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void array<X>::setRandom(X lower, X upper)
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

template<class X>
inline
void array<X>::setArange()
{
  std::iota(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void array<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void array<X>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void array<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// =================================================================================================
// copy to target
// =================================================================================================

template<class X>
template<class Iterator>
inline
void array<X>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline
void array<X>::copyTo(Iterator first, Iterator last) const
{
  assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// =================================================================================================
// bound check
// =================================================================================================

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a, T b) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;
  if ( b >= static_cast<T>(mShape[1]) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a, T b, T c) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;
  if ( b >= static_cast<T>(mShape[1]) ) return false;
  if ( c >= static_cast<T>(mShape[2]) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
    if ( c < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a, T b, T c, T d) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;
  if ( b >= static_cast<T>(mShape[1]) ) return false;
  if ( c >= static_cast<T>(mShape[2]) ) return false;
  if ( d >= static_cast<T>(mShape[3]) ) return false;

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

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a, T b, T c, T d, T e) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;
  if ( b >= static_cast<T>(mShape[1]) ) return false;
  if ( c >= static_cast<T>(mShape[2]) ) return false;
  if ( d >= static_cast<T>(mShape[3]) ) return false;
  if ( e >= static_cast<T>(mShape[4]) ) return false;

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

template<class X>
template<typename T>
inline
bool array<X>::inBounds(T a, T b, T c, T d, T e, T f) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(mShape[0]) ) return false;
  if ( b >= static_cast<T>(mShape[1]) ) return false;
  if ( c >= static_cast<T>(mShape[2]) ) return false;
  if ( d >= static_cast<T>(mShape[3]) ) return false;
  if ( e >= static_cast<T>(mShape[4]) ) return false;
  if ( f >= static_cast<T>(mShape[5]) ) return false;

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

template<class X>
inline
array<X> array<X>::operator- () const
{
  array<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::operator+ () const
{
  array<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = mData[i];

  return out;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline
array<X>& array<X>::operator*= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator/= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator+= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator-= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator*= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator/= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator+= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X>& array<X>::operator-= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<class X>
inline
array<X> array<X>::abs() const
{
  array<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = std::abs(mData[i]);

  return out;
}

// =================================================================================================
// norm
// =================================================================================================

template<class X>
inline
X array<X>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// return the indices that would sort the array
// =================================================================================================

template<class X>
inline
array<size_t> array<X>::argsort(bool ascending) const
{
  return array<size_t>::Copy(shape(), cppmat::argsort(mData, ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<class X>
inline
size_t array<X>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
size_t array<X>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<class X>
inline
X array<X>::min() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::min(size_t axis) const
{
  // check input
  assert( axis < mRank );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(shape(),axis), max());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::min(out[ni], mData[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::min(int axis) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( axis  <      static_cast<int>(mRank) );
  assert( axis >= -1 * static_cast<int>(mRank) );

  // get number of dimensions as integer
  int n = static_cast<int>(mRank);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return min(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::min(const std::vector<int> &axes_in) const
{
  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = cppmat::Private::sort_axes(axes_in, static_cast<int>(mRank), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.min(axis);

  return out;
}

// =================================================================================================
// maximum
// =================================================================================================

template<class X>
inline
X array<X>::max() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::max(size_t axis) const
{
  // check input
  assert( axis < mRank );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(shape(),axis), min());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::max(out[ni], mData[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::max(int axis) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( axis  <      static_cast<int>(mRank) );
  assert( axis >= -1 * static_cast<int>(mRank) );

  // get number of dimensions as integer
  int n = static_cast<int>(mRank);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return max(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::max(const std::vector<int> &axes_in) const
{
  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = cppmat::Private::sort_axes(axes_in, static_cast<int>(mRank), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.max(axis);

  return out;
}

// =================================================================================================
// sum
// =================================================================================================

template<class X>
inline
X array<X>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::sum(size_t axis) const
{
  // zero-initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Zero(del(shape(),axis));

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] += mData[i];
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::sum(int axis) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( axis  <      static_cast<int>(mRank) );
  assert( axis >= -1 * static_cast<int>(mRank) );

  // get number of dimensions as integer
  int n = static_cast<int>(mRank);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return sum(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::sum(const std::vector<int> &axes_in) const
{
  // check rank
  assert( axes_in.size() < mRank );

  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = cppmat::Private::sort_axes(axes_in, static_cast<int>(mRank), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.sum(axis);

  return out;
}

// =================================================================================================
// mean
// =================================================================================================

template<class X>
inline
double array<X>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::mean(size_t axis) const
{
  array<X> weights = array<X>::Ones(shape());

  return (*this).sum(axis) / weights.sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::mean(int axis) const
{
  array<X> weights = array<X>::Ones(shape());

  return (*this).sum(axis) / weights.sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::mean(const std::vector<int> &axes) const
{
  array<X> weights = array<X>::Ones(shape());

  return (*this).sum(axes) / weights.sum(axes);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<class X>
inline
double array<X>::average(const array<X> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::average(const array<X> &weights, size_t axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::average(const array<X> &weights, int axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> array<X>::average(
  const array<X> &weights, const std::vector<int> &axes, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axes) / weights.sum(axes);
  else        return (weights*(*this)).sum(axes);
}

// =================================================================================================
// return array of booleans, based on condition
// =================================================================================================

template<class X>
inline
array<int> array<X>::equal(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::not_equal(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::greater(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::greater_equal(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::less(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::less_equal(const X &D) const
{
  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::equal(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::not_equal(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::greater(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::greater_equal(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::less(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<int> array<X>::less_equal(const array<X> &D) const
{
  assert( shape() == D.shape() );
  assert( size () == D.size () );

  array<int> out = array<int>::Zero(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D[i] )
      out[i] = 1;

  return out;
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<class X>
inline
std::vector<size_t> array<X>::where() const
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

template<class X>
inline
size_t array<X>::where(size_t index) const
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

template<class X>
inline
size_t array<X>::where(int index) const
{
  int nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  assert( index < nnz && index >= -nnz );

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

template<class X>
inline
std::ostream& operator<<(std::ostream& out, const array<X>& src)
{
  if ( src.size() == 0 ) return out;

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

template<class X>
inline
bool operator!= (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
bool operator== (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<class X>
inline
array<X> operator* (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator/ (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator+ (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator- (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator* (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator/ (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator+ (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator- (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator* (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator/ (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator+ (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
array<X> operator- (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// minimum/maximum from two arrays of equal shape
// =================================================================================================

template<class X> array<X> min(const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = std::min(A[i], B[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> array<X> max(const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.rank () == B.rank () );
  assert( A.size () == B.size () );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = std::max(A[i], B[i]);

  return C;
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

