/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_DIAGONAL_MATRIX_HPP
#define CPPMAT_VAR_DIAGONAL_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace diagonal {

// =================================================================================================
// constructors
// =================================================================================================

template<typename X>
inline
matrix<X>::matrix()
{
  mZero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X>::matrix(size_t m, size_t n)
{
  mZero[0] = static_cast<X>(0);

  Assert( m == n );

  resize(m,n);
}

// =================================================================================================
// constructors: copy from own class (with different type)
// =================================================================================================

template<typename X>
template<typename U, typename V>
inline
matrix<X>::matrix(const matrix<U> &A)
{
  mZero[0] = static_cast<X>(0);

  resize(A.shape(0), A.shape(1));

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = static_cast<X>(A[i]);
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<typename X>
template<size_t m, size_t n>
inline
matrix<X>::matrix(const cppmat::tiny::diagonal::matrix<X,m,n> &A) : cppmat::diagonal::matrix<X>(m,n)
{
  setCopy(A.begin(), A.end());
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<typename X>
template<size_t m, size_t n>
inline
matrix<X>::matrix(const cppmat::view::diagonal::matrix<X,m,n> &A) : cppmat::diagonal::matrix<X>(m,n)
{
  setCopy(A.begin(), A.end());
}

// =================================================================================================
// named constructors
// =================================================================================================

template<typename X>
inline
matrix<X> matrix<X>::Random(size_t m, size_t n, X lower, X upper)
{
  matrix<X> out(m,n);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::Arange(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::Zero(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::Ones(size_t m, size_t n)
{
  matrix<X> out(m,n);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::Constant(size_t m, size_t n, X D)
{
  matrix<X> out(m,n);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::Copy(size_t m, size_t n, const std::vector<X> &D)
{
  matrix<X> out(m,n);

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
matrix<X> matrix<X>::Copy(size_t m, size_t n, Iterator first)
{
  matrix<X> out(m,n);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
matrix<X> matrix<X>::Copy(size_t m, size_t n, Iterator first, Iterator last)
{
  matrix<X> out(m,n);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// named constructor: copy from dense storage
// =================================================================================================

template<typename X>
template<typename Iterator>
inline
matrix<X> matrix<X>::CopyDense(size_t m, size_t n, Iterator first)
{
  Assert( m == n );

  matrix<X> out(m,n);

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
matrix<X> matrix<X>::CopyDense(size_t m, size_t n, Iterator first, Iterator last)
{
  Assert( m == n );

  matrix<X> out(m,n);

  out.setCopyDense(first,last);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<typename X>
template<typename U, typename V>
inline
matrix<X>::operator std::vector<U> () const
{
  std::vector<U> out(mSize);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = static_cast<U>(mData[i]);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<typename X>
inline
void matrix<X>::resize(size_t m, size_t n)
{
  Assert( m == n );

  // store old size
  size_t size = mSize;

  // copy to class member
  N = m;

  // set number of dimensions and total size
  mSize = N;

  // resize data container
  if ( mSize != size ) mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
void matrix<X>::resize(size_t m, size_t n, const X &D)
{
  Assert( m == n );

  // store old size
  size_t size = mSize;

  // copy to class member
  N = m;

  // set number of dimensions and total size
  mSize = N;

  // resize data container
  if ( mSize != size ) mData.resize(mSize, D);
}

// =================================================================================================
// modify bounds check
// =================================================================================================

template<typename X>
inline
void matrix<X>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<typename X>
inline
size_t matrix<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
size_t matrix<X>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
size_t matrix<X>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  Assert( i  <      static_cast<int>(mRank) );
  Assert( i >= -1 * static_cast<int>(mRank) );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
size_t matrix<X>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  Assert( i < mRank );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
std::vector<size_t> matrix<X>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::fill(shape.begin(), shape.end(), N);

  return shape;
}

// =================================================================================================
// get dimensions using a different return type
// =================================================================================================

template<typename X>
template<typename U>
inline
U matrix<X>::size() const
{
  return static_cast<U>(size());
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename U>
inline
U matrix<X>::rank() const
{
  return static_cast<U>(rank());
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename U>
inline
U matrix<X>::shape(int i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename U>
inline
U matrix<X>::shape(size_t i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename U>
inline
std::vector<U> matrix<X>::shape() const
{
  std::vector<size_t> A = shape();

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<typename X>
inline
X& matrix<X>::operator[](size_t i)
{
  Assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
const X& matrix<X>::operator[](size_t i) const
{
  Assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<typename X>
inline
X& matrix<X>::operator()(int a, int b)
{
  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );
  Assert( ( b < n && b >= -n ) or mPeriodic );

  if ( a != b ) return mZero[0];

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return mData[A];
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
const X& matrix<X>::operator()(int a, int b) const
{
  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );
  Assert( ( b < n && b >= -n ) or mPeriodic );

  if ( a != b ) return mZero[0];

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return mData[A];
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T, typename S>
inline
X& matrix<X>::operator()(T a, T b)
{
  Assert( a < N );
  Assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T, typename S>
inline
const X& matrix<X>::operator()(T a, T b) const
{
  Assert( a < N );
  Assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<typename X>
inline
size_t matrix<X>::compress(int a, int b) const
{
  Assert( a == b );

  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return A;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T, typename S>
inline
size_t matrix<X>::compress(T a, T b) const
{
  Assert( a < N  );
  Assert( b < N  );
  Assert( a == b );

  return a;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<typename X>
inline
std::vector<size_t> matrix<X>::decompress(size_t i) const
{
  // check input
  Assert( i < mSize );

  // allocate matrix-index
  std::vector<size_t> idx(mRank);

  // reconstruct
  std::fill(idx.begin(), idx.end(), i);

  return idx;
}

// =================================================================================================
// midpoint
// =================================================================================================

template<typename X>
inline
std::vector<size_t> matrix<X>::midpoint() const
{
  // get shape
  std::vector<size_t> mid = shape();

  // check odd-sized
  for ( auto &i : mid )
    if ( i%2 == 0 )
      throw std::domain_error("cppmat::matrix<X>::midpoint: Must be odd shaped");

  // midpoint
  for ( auto &i : mid )
    i = (i-1)/2;

  return mid;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
size_t matrix<X>::midpoint(size_t axis) const
{
  // get shape
  size_t mid = shape(axis);

  // check odd-sized
  if ( mid%2 == 0 )
    throw std::domain_error("cppmat::matrix<X>::midpoint: Must be odd shaped");

  // midpoint
  mid = (mid-1)/2;

  return mid;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<typename X>
inline
X* matrix<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
const X* matrix<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<typename X>
inline
auto matrix<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
auto matrix<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
auto matrix<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
auto matrix<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<typename X>
inline
auto matrix<X>::index(size_t i)
{
  Assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
auto matrix<X>::index(size_t i) const
{
  Assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<typename X>
inline
auto matrix<X>::item(int a, int b)
{
  Assert( a == b );

  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return this->begin() + A;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
auto matrix<X>::item(int a, int b) const
{
  Assert( a == b );

  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return this->begin() + A;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T, typename S>
inline
auto matrix<X>::item(T a, T b)
{
  Assert( a < N  );
  Assert( b < N  );
  Assert( a == b );

  return this->begin() + a;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T, typename S>
inline
auto matrix<X>::item(T a, T b) const
{
  Assert( a < N  );
  Assert( b < N  );
  Assert( a == b );

  return this->begin() + a;
}

// =================================================================================================
// initialize
// =================================================================================================

template<typename X>
inline
void matrix<X>::setRandom(X lower, X upper)
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

template<typename X>
inline
void matrix<X>::setArange()
{
  std::iota(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
void matrix<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
void matrix<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
void matrix<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<class Iterator>
inline
void matrix<X>::setCopy(Iterator first)
{
  std::copy(first, first+mSize, begin());
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<class Iterator>
inline
void matrix<X>::setCopy(Iterator first, Iterator last)
{
  Assert( mSize == static_cast<size_t>(last-first) );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
void matrix<X>::setCopyDense(Iterator first)
{
  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        Assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
void matrix<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  Assert( N*N == static_cast<size_t>(last-first) );

  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      if ( i !=j )
        Assert( !first[i*N+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < N ; ++i )
    mData[i] = first[i*N+i];
}

// =================================================================================================
// copy to target
// =================================================================================================

template<typename X>
template<class Iterator>
inline
void matrix<X>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<class Iterator>
inline
void matrix<X>::copyTo(Iterator first, Iterator last) const
{
  Assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first) const
{
  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<class Iterator>
inline
void matrix<X>::copyToDense(Iterator first, Iterator last) const
{
  Assert( N*N == static_cast<size_t>(last-first) );

  UNUSED(last);

  for ( size_t i = 0 ; i < N ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
	  if ( i == j ) first[i*N+j] = mData[i];
	  else          first[i*N+j] = static_cast<X>(0);
	}
  }
}

// =================================================================================================
// bound check
// =================================================================================================

template<typename X>
template<typename T>
inline
bool matrix<X>::inBounds(T a) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(N) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename T>
inline
bool matrix<X>::inBounds(T a, T b) const
{
  if ( mPeriodic ) return true;

  if ( a >= static_cast<T>(N) ) return false;
  if ( b >= static_cast<T>(N) ) return false;

  if ( std::numeric_limits<T>::is_signed )
  {
    if ( a < 0 ) return false;
    if ( b < 0 ) return false;
  }

  return true;
}

// =================================================================================================
// sign change
// =================================================================================================

template<typename X>
inline
matrix<X> matrix<X>::operator- () const
{
  matrix<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = -mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> matrix<X>::operator+ () const
{
  matrix<X> out(shape());

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = mData[i];

  return out;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<typename X>
inline
matrix<X>& matrix<X>::operator*= (const matrix<X> &B)
{
  Assert( shape() == B.shape() );
  Assert( rank () == B.rank () );
  Assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X>& matrix<X>::operator+= (const matrix<X> &B)
{
  Assert( shape() == B.shape() );
  Assert( rank () == B.rank () );
  Assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X>& matrix<X>::operator-= (const matrix<X> &B)
{
  Assert( shape() == B.shape() );
  Assert( rank () == B.rank () );
  Assert( size () == B.size () );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X>& matrix<X>::operator*= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X>& matrix<X>::operator/= (X B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// =================================================================================================
// absolute value
// =================================================================================================

template<typename X>
inline
matrix<X> matrix<X>::abs() const
{
  matrix<X> out(N, N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = std::abs(mData[i]);

  return out;
}

// =================================================================================================
// norm
// =================================================================================================

template<typename X>
inline
X matrix<X>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// return the indices that would sort the matrix
// =================================================================================================

template<typename X>
inline
matrix<size_t> matrix<X>::argsort(bool ascending) const
{
  return matrix<size_t>::Copy(N, N, cppmat::argsort(mData, ascending));
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<typename X>
inline
size_t matrix<X>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
size_t matrix<X>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<typename X>
inline
X matrix<X>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<typename X>
inline
X matrix<X>::max() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<typename X>
inline
X matrix<X>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// =================================================================================================
// mean
// =================================================================================================

template<typename X>
inline
double matrix<X>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(N*N);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<typename X>
inline
double matrix<X>::average(const matrix<X> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// return array of booleans, based on condition
// =================================================================================================

template<typename X>
inline
matrix<int> matrix<X>::equal(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::not_equal(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::greater(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::greater_equal(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::less(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::less_equal(const X &D) const
{
  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::equal(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] == D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::not_equal(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::greater(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] > D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::greater_equal(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] >= D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::less(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] < D[i] )
      out[i] = 1;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<int> matrix<X>::less_equal(const matrix<X> &D) const
{
  Assert( shape() == D.shape() );
  Assert( size () == D.size () );

  matrix<int> out = matrix<int>::Zero(N,N);

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] <= D[i] )
      out[i] = 1;

  return out;
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<typename X>
inline
std::vector<size_t> matrix<X>::where() const
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

template<typename X>
inline
size_t matrix<X>::where(size_t index) const
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

template<typename X>
inline
size_t matrix<X>::where(int index) const
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

template<typename X>
inline
std::ostream& operator<<(std::ostream& out, const matrix<X>& src)
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
// equality operators
// =================================================================================================

template<typename X>
inline
bool operator!= (const matrix<X> &A, const matrix<X> &B)
{
  Assert( A.shape() == B.shape() );
  Assert( A.rank () == B.rank () );
  Assert( A.size () == B.size () );

  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
bool operator== (const matrix<X> &A, const matrix<X> &B)
{
  Assert( A.shape() == B.shape() );
  Assert( A.rank () == B.rank () );
  Assert( A.size () == B.size () );

  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================
// arithmetic operators: external
// =================================================================================================

template<typename X>
inline
matrix<X> operator* (const matrix<X> &A, const matrix<X> &B)
{
  Assert( A.shape() == B.shape() );
  Assert( A.rank () == B.rank () );
  Assert( A.size () == B.size () );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B)
{
  Assert( A.shape() == B.shape() );
  Assert( A.rank () == B.rank () );
  Assert( A.size () == B.size () );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> operator- (const matrix<X> &A, const matrix<X> &B)
{
  Assert( A.shape() == B.shape() );
  Assert( A.rank () == B.rank () );
  Assert( A.size () == B.size () );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> operator* (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> operator/ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
matrix<X> operator* (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

