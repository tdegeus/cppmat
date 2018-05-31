/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_SYMMETRIC_MATRIX_CPP
#define CPPMAT_MAP_SYMMETRIC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace symmetric {

// =================================================================================================
// return size without constructing
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::Size()
{
  return (N+1)*N/2;
}

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix()
{
}

// =================================================================================================
// constructors: map external pointer
// =================================================================================================

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::matrix(const X *A)
{
  mData = A;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X, size_t M, size_t N>
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

template<class X, size_t M, size_t N>
inline
matrix<X,M,N>::operator std::vector<X> () const
{
  std::vector<X> out(mSize);

  std::copy(begin(), end(), out.begin());

  return out;
}

// =================================================================================================
// modify bounds check
// =================================================================================================

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setPeriodic(bool periodic)
{
  mPeriodic = periodic;
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
const X& matrix<X,M,N>::operator()(int a, int b) const
{
  int n = static_cast<int>(N);

  assert( ( a < n && a >= -n ) or mPeriodic );
  assert( ( b < n && b >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );
  size_t B = static_cast<size_t>( (n+(b%n)) % n );

  if (A <= B) return mData[ A*N - (A-1)*A/2 + B - A ];
  else        return mData[ B*N - (B-1)*B/2 + A - B ];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
const X& matrix<X,M,N>::operator()(T a, T b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return mData[ a*N - (a-1)*a/2 + b - a ];
  else        return mData[ b*N - (b-1)*b/2 + a - b ];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline
size_t matrix<X,M,N>::compress(int a, int b) const
{
  int n = static_cast<int>(N);

  assert( ( a < n && a >= -n ) or mPeriodic );
  assert( ( b < n && b >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );
  size_t B = static_cast<size_t>( (n+(b%n)) % n );

  if (A <= B) return A*N - (A-1)*A/2 + B - A;
  else        return B*N - (B-1)*B/2 + A - B;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
size_t matrix<X,M,N>::compress(T a, T b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return a*N - (a-1)*a/2 + b - a;
  else        return b*N - (b-1)*b/2 + a - b;
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
  idx[0] = 0;
  size_t keyafter;
  do
  {
    idx[0]++;
    keyafter = idx[0] * N - (idx[0] - 1) * idx[0] / 2;
  } while ( i >= keyafter );
  idx[0]--;
  idx[1] = N - keyafter + i;

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t M, size_t N>
inline
const X* matrix<X,M,N>::data() const
{
  return mData;
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::begin() const
{
  return mData;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline
auto matrix<X,M,N>::end() const
{
  return mData + mSize;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

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
auto matrix<X,M,N>::item(int a, int b) const
{
  int n = static_cast<int>(N);

  assert( ( a < n && a >= -n ) or mPeriodic );
  assert( ( b < n && b >= -n ) or mPeriodic );

  size_t A = static_cast<size_t>( (n+(a%n)) % n );
  size_t B = static_cast<size_t>( (n+(b%n)) % n );

  if (A <= B) return begin() + ( A*N - (A-1)*A/2 + B - A );
  else        return begin() + ( B*N - (B-1)*B/2 + A - B );
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<typename T, typename S>
inline
auto matrix<X,M,N>::item(T a, T b) const
{
  assert( a < N );
  assert( b < N );

  if (a <= b) return begin() + ( a*N - (a-1)*a/2 + b - a );
  else        return begin() + ( b*N - (b-1)*b/2 + a - b );
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X, size_t M, size_t N>
inline
void matrix<X,M,N>::setMap(const X *D)
{
  mData = D;
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
  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      first[i*N+j] = first[j*N+i] = mData[i*N-(i-1)*i/2+j-i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
template<class Iterator>
inline
void matrix<X,M,N>::copyToDense(Iterator first, Iterator last) const
{
  assert( N*N == static_cast<size_t>(last-first) );

  UNUSED(last);

  for ( size_t i = 0 ; i < N ; ++i )
    for ( size_t j = i ; j < N ; ++j )
      first[i*N+j] = first[j*N+i] = mData[i*N-(i-1)*i/2+j-i];
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
  X out = static_cast<X>(0);

  for ( size_t a = 0 ; a < N ; ++a ) {
    for ( size_t b = a ; b < N ; ++b ) {
      if (a == b) out += mData[ a*N - (a-1)*a/2 + b - a ];
      else        out += mData[ a*N - (a-1)*a/2 + b - a ] * static_cast<X>(2);
    }
  }

  return out;
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

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

