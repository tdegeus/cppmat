/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_DIAGONAL_MATRIX_HPP
#define CPPMAT_FIX_DIAGONAL_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace diagonal {

// =================================================================================================
// return size without constructing
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::Size()
{
  return N;
}

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix()
{
  mZero[0] = static_cast<X>(0);
}

// =================================================================================================
// constructors: copy from own class (with different type)
// =================================================================================================

template<class X, size_t M, size_t N>
template<typename U, typename V>
inline
matrix<X,M,N>::matrix(const matrix<U,M,N> &A)
{
  mZero[0] = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = static_cast<X>(A[i]);
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::diagonal::matrix<X> &A) : cppmat::tiny::diagonal::matrix<X,M,N>()
{
  assert( N == A.shape(0) );
  assert( N == A.shape(1) );

  setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const cppmat::view::diagonal::matrix<X,M,N> &A) : cppmat::tiny::diagonal::matrix<X,M,N>()
{
  setCopy(A.begin(), A.end());
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Random(X lower, X upper)
{
  matrix<X,M,N> out;

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Arange()
{
  matrix<X,M,N> out;

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Zero()
{
  matrix<X,M,N> out;

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Ones()
{
  matrix<X,M,N> out;

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Constant(X D)
{
  matrix<X,M,N> out;

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Copy(const std::vector<X> &D)
{
  matrix<X,M,N> out;

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
matrix<X,M,N> matrix<X,M,N>::Copy(Iterator first)
{
  matrix<X,M,N> out;

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
matrix<X,M,N> matrix<X,M,N>::Copy(Iterator first, Iterator last)
{
  matrix<X,M,N> out;

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// named constructor: copy from dense storage
// =================================================================================================

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
matrix<X,M,N> matrix<X,M,N>::CopyDense(Iterator first)
{
  matrix<X,M,N> out;

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
matrix<X,M,N> matrix<X,M,N>::CopyDense(Iterator first, Iterator last)
{
  matrix<X,M,N> out;

  out.setCopyDense(first,last);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::operator std::vector<X> () const
{
  std::vector<X> out(mSize);

  std::copy(begin(), end(), out.begin());

  return out;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  assert( i  <      static_cast<int>(mRank) );
  assert( i >= -1 * static_cast<int>(mRank) );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  assert( i < mRank );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::fill(shape.begin(), shape.end(), N);

  return shape;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X& matrix<X,M,N>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
const X& matrix<X,M,N>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X& matrix<X,M,N>::operator()(int a, int b)
{
  int n = static_cast<int>(N);

  if ( a != b ) return mZero[0];

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return mData[A];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
const X& matrix<X,M,N>::operator()(int a, int b) const
{
  int n = static_cast<int>(N);

  if ( a != b ) return mZero[0];

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return mData[A];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
X& matrix<X,M,N>::operator()(T a, T b)
{
  assert( a < N );
  assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
const X& matrix<X,M,N>::operator()(T a, T b) const
{
  assert( a < N );
  assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::compress(int a, int b) const
{
  assert( a == b );

  int n = static_cast<int>(N);

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
size_t matrix<X,M,N>::compress(T a, T b) const
{
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return a;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate matrix-index
  std::vector<size_t> idx(mRank);

  // reconstruct
  std::fill(idx.begin(), idx.end(), i);

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X* matrix<X,M,N>::data()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
const X* matrix<X,M,N>::data() const
{
  return std::begin(mData);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::begin()
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::begin() const
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::end()
{
  return std::begin(mData) + mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::end() const
{
  return std::begin(mData) + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::item(int a, int b)
{
  assert( a == b );

  int n = static_cast<int>(N);

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return this->begin() + A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::item(int a, int b) const
{
  assert( a == b );

  int n = static_cast<int>(N);

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return this->begin() + A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
auto matrix<X,M,N>::item(T a, T b)
{
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return this->begin() + a;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
auto matrix<X,M,N>::item(T a, T b) const
{
  assert( a < N  );
  assert( b < N  );
  assert( a == b );

  return this->begin() + a;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setRandom(X lower, X upper)
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

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setArange()
{
  std::iota(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
void matrix<X,M,N>::setCopyDense(Iterator first)
{
  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename Iterator>
inline
void matrix<X,M,N>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( N*N == static_cast<size_t>(last-first) );

  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
}

// =================================================================================================
// copy to target
// =================================================================================================

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyTo(Iterator first, Iterator last) const
{
  assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyToDense(Iterator first) const
{
  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyToDense(Iterator first, Iterator last) const
{
  assert( N*N == static_cast<size_t>(last-first) );

  UNUSED(last);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
}

// =================================================================================================
// sign change
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::operator- () const
{
  matrix<X,M,N> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::operator+ () const
{
  matrix<X,M,N> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = mData[i];

  return out;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>& matrix<X,M,N>::operator*= (const matrix<X,M,N> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>& matrix<X,M,N>::operator+= (const matrix<X,M,N> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>& matrix<X,M,N>::operator-= (const matrix<X,M,N> &B)
{
  assert( shape() == B.shape() );
  assert( rank () == B.rank () );
  assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>& matrix<X,M,N>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>& matrix<X,M,N>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::abs() const
{
  matrix<X,M,N> out;

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = std::abs(mData[i]);

  return out;
}

// =================================================================================================
// norm
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X matrix<X,M,N>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// return the indices that would sort the matrix
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<size_t,M,N> matrix<X,M,N>::argsort(bool ascending) const
{
  return matrix<size_t,M,N>::Copy(cppmat::argsort(mData, ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X matrix<X,M,N>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X matrix<X,M,N>::max() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<class X, size_t M, size_t N>
inline
X matrix<X,M,N>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// =================================================================================================
// mean
// =================================================================================================

template<class X, size_t M, size_t N>
inline
double matrix<X,M,N>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(N*N);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<class X, size_t M, size_t N>
inline
double matrix<X,M,N>::average(const matrix<X,M,N> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// return array of booleans, based on condition
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::equal(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::not_equal(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::greater(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::greater_equal(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::less(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::less_equal(const X &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::equal(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::not_equal(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::greater(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::greater_equal(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::less(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<int,M,N> matrix<X,M,N>::less_equal(const matrix<X,M,N> &D) const
{
  matrix<int,M,N> out = matrix<int,M,N>::Zero();

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D[i] )
      out[i] = 1;

  return out;
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<class X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::where() const
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

template<class X, size_t M, size_t N>
inline
std::ostream& operator<<(std::ostream& out, const matrix<X,M,N>& src)
{
  auto w = out.width();
  auto p = out.precision();

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

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator* (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator- (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator* (const matrix<X,M,N> &A, const X &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const X &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
matrix<X,M,N> operator* (const X &A, const matrix<X,M,N> &B)
{
  matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

