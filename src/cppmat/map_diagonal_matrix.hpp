/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_DIAGONAL_MATRIX_HPP
#define CPPMAT_MAP_DIAGONAL_MATRIX_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace diagonal {

// =================================================================================================
// return size without constructing
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::Size()
{
  return N;
}

// =================================================================================================
// constructors
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix()
{
  mZero[0] = static_cast<X>(0);
}

// =================================================================================================
// constructors: map external pointer
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const X *A)
{
  mZero[0] = static_cast<X>(0);

  mData = A;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
matrix<X,M,N> matrix<X,M,N>::Map(const X *D)
{
  matrix<X,M,N> out;

  out.setMap(D);

  return out;
}

// =================================================================================================
// return plain storage as vector
// =================================================================================================

template<typename X, size_t M, size_t N>
template<typename U, typename V>
inline
matrix<X,M,N>::operator std::vector<U> () const
{
  std::vector<U> out(mSize);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out[i] = static_cast<U>(mData[i]);

  return out;
}

// =================================================================================================
// modify bounds check
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
void matrix<X,M,N>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::rank() const
{
  return mRank;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::shape(int i) const
{
  // check axis: (0,1,...,rank-1) or (-1,-2,...,-rank)
  Assert( i  <      static_cast<int>(mRank) );
  Assert( i >= -1 * static_cast<int>(mRank) );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::shape(size_t i) const
{
  // check axis: (0,1,...,rank-1)
  Assert( i < mRank );

  // return shape
  return N;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::shape() const
{
  std::vector<size_t> shape(mRank);

  std::fill(shape.begin(), shape.end(), N);

  return shape;
}

// =================================================================================================
// get dimensions using a different return type
// =================================================================================================

template<typename X, size_t M, size_t N>
template<typename U>
inline
U matrix<X,M,N>::size() const
{
  return static_cast<U>(size());
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename U>
inline
U matrix<X,M,N>::rank() const
{
  return static_cast<U>(rank());
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename U>
inline
U matrix<X,M,N>::shape(int i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename U>
inline
U matrix<X,M,N>::shape(size_t i) const
{
  return static_cast<U>(shape(i));
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename U>
inline
std::vector<U> matrix<X,M,N>::shape() const
{
  std::vector<size_t> A = shape();

  std::vector<U> B(A.size());

  std::copy(A.begin(), A.end(), B.begin());

  return B;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
const X& matrix<X,M,N>::operator[](size_t i) const
{
  Assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
const X& matrix<X,M,N>::operator()(int a, int b) const
{
  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );
  Assert( ( b < n && b >= -n ) or mPeriodic );

  if ( a != b ) return mZero[0];

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return mData[A];
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename T, typename S>
inline
const X& matrix<X,M,N>::operator()(T a, T b) const
{
  Assert( a < N );
  Assert( b < N );

  if (a == b) return mData[a];
  else        return mZero[0];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::compress(int a, int b) const
{
  Assert( a == b );

  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return A;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename T, typename S>
inline
size_t matrix<X,M,N>::compress(T a, T b) const
{
  Assert( a < N  );
  Assert( b < N  );
  Assert( a == b );

  return a;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::decompress(size_t i) const
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

template<typename X, size_t M, size_t N>
inline
std::vector<size_t> matrix<X,M,N>::midpoint() const
{
  // get shape
  std::vector<size_t> mid = shape();

  // check odd-sized
  for ( auto &i : mid )
    if ( i%2 == 0 )
      throw std::domain_error("cppmat::matrix<X,M,N>::midpoint: Must be odd shaped");

  // midpoint
  for ( auto &i : mid )
    i = (i-1)/2;

  return mid;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::midpoint(size_t axis) const
{
  // get shape
  size_t mid = shape(axis);

  // check odd-sized
  if ( mid%2 == 0 )
    throw std::domain_error("cppmat::matrix<X,M,N>::midpoint: Must be odd shaped");

  // midpoint
  mid = (mid-1)/2;

  return mid;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
const X* matrix<X,M,N>::data() const
{
  return mData;
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
auto matrix<X,M,N>::begin() const
{
  return mData;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
auto matrix<X,M,N>::end() const
{
  return mData + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
auto matrix<X,M,N>::index(size_t i) const
{
  Assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
auto matrix<X,M,N>::item(int a, int b) const
{
  Assert( a == b );

  int n = static_cast<int>(N);

  Assert( ( a < n && a >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );

  return this->begin() + A;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<typename T, typename S>
inline
auto matrix<X,M,N>::item(T a, T b) const
{
  Assert( a < N  );
  Assert( b < N  );
  Assert( a == b );

  return this->begin() + a;
}

// =================================================================================================
// initialize
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
void matrix<X,M,N>::setMap(const X *D)
{
  mData = D;
}

// =================================================================================================
// copy to target
// =================================================================================================

template<typename X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyTo(Iterator first) const
{
  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyTo(Iterator first, Iterator last) const
{
  Assert( mSize == static_cast<size_t>(last-first) );

  UNUSED(last);

  std::copy(begin(), end(), first);
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
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

template<typename X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyToDense(Iterator first, Iterator last) const
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

template<typename X, size_t M, size_t N>
template<typename T>
inline
bool matrix<X,M,N>::inBounds(T a) const
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

template<typename X, size_t M, size_t N>
template<typename T>
inline
bool matrix<X,M,N>::inBounds(T a, T b) const
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
// norm
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
X matrix<X,M,N>::norm() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += std::abs(mData[i]);

  return out;
}

// =================================================================================================
// location of the minimum/maximum
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::argmin() const
{
  return std::min_element(begin(), end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::argmax() const
{
  return std::max_element(begin(), end()) - begin();
}

// =================================================================================================
// minimum
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
X matrix<X,M,N>::min() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// maximum
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
X matrix<X,M,N>::max() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// sum
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
X matrix<X,M,N>::sum() const
{
  return std::accumulate(begin(), end(), static_cast<X>(0));
}

// =================================================================================================
// mean
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
double matrix<X,M,N>::mean() const
{
  return static_cast<double>(sum())/static_cast<double>(N*N);
}

// =================================================================================================
// weighted average
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
double matrix<X,M,N>::average(const matrix<X,M,N> &weights, bool norm) const
{
  if ( norm ) return static_cast<double>((weights*(*this)).sum())/static_cast<double>(weights.sum());
  else        return static_cast<double>((weights*(*this)).sum());
}

// =================================================================================================
// find the plain storage indices of all non-zero entries
// =================================================================================================

template<typename X, size_t M, size_t N>
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

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::where(size_t index) const
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

template<typename X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::where(int index) const
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

template<typename X, size_t M, size_t N>
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
// equality operators
// =================================================================================================

template<typename X, size_t M, size_t N>
inline
bool operator!= (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return true;

  return false;
}

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
inline
bool operator== (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  for ( size_t i = 0 ; i < A.size() ; ++i )
    if ( A[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

