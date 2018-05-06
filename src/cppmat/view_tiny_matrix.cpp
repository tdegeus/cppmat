/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TINY_MATRIX_CPP
#define CPPMAT_VIEW_TINY_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "view_tiny_matrix.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t M, size_t N>
inline matrix<X,M,N> matrix<X,M,N>::Map(const X *D)
{
  // call basic constructor
  matrix<X,M,N> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::Size()
{
  return M*N;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t M, size_t N>
inline void matrix<X,M,N>::setMap(const X *D)
{
  mData = D;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::ndim() const
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::rows() const
{
  return M;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::cols() const
{
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::shape(int i) const
{
  // check axis: (0,1) or (-1,-2)
  assert( i  <  2 );
  assert( i >= -2 );

  // correct periodic index
  i = ( 2 + (i%2) ) % 2;

  // return shape
  if ( i == 0 ) return M;
  else          return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::shape(size_t i) const
{
  // check axis: (0,1)
  assert( i < 2 );

  // return shape
  if ( i == 0 ) return M;
  else          return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline std::vector<size_t> matrix<X,M,N>::shape() const
{
  std::vector<size_t> ret(2);

  ret[0] = M;
  ret[1] = N;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline std::vector<size_t> matrix<X,M,N>::strides(bool bytes) const
{
  std::vector<size_t> ret(2);

  ret[0] = N;
  ret[1] = 1;

  if ( bytes ) {
    ret[0] *= sizeof(X);
    ret[1] *= sizeof(X);
  }

  return ret;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X, size_t M, size_t N>
inline const X& matrix<X,M,N>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline const X& matrix<X,M,N>::operator()(size_t a) const
{
  assert( a < M );

  return mData[a*N];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline const X& matrix<X,M,N>::operator()(size_t a, size_t b) const
{
  assert( a < M );
  assert( b < N );

  return mData[a*N+b];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X, size_t M, size_t N>
template<class Iterator>
inline const X& matrix<X,M,N>::at(Iterator first, Iterator last) const
{
  // check input
  assert( last-first  > 0 );
  assert( last-first <= 2 );

  // suppress compiler warning
  UNUSED(last);

  // index
  int a = first[0];
  int b = 0;

  // optional index
  if ( last-first == 2 ) b = first[1];

  return mData[a*N+b];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::compress(size_t a) const
{
  assert( a < M );

  return a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::compress(size_t a, size_t b) const
{
  assert( a < M );
  assert( b < N );

  return a*N+b;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X, size_t M, size_t N>
inline std::vector<size_t> matrix<X,M,N>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(2);

  // reconstruct
  idx[1] = i % N;
  idx[0] = ( i - idx[1] ) / N;

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t M, size_t N>
inline const X* matrix<X,M,N>::data() const
{
  return std::begin(mData);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::begin() const
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::end() const
{
  return std::begin(mData) + mSize;
}

// =================================================================================================
// iterators : beginRow() and endRow()
// =================================================================================================

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::beginRow(size_t a) const
{
  assert( a < M );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::endRow(size_t a) const
{
  assert( a < M );

  return begin() + (a+1)*N;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::item(size_t a) const
{
  assert( a < M );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline auto matrix<X,M,N>::item(size_t a, size_t b) const
{
  assert( a < M );
  assert( b < N );

  return begin() + a*N+b;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const matrix<X,M,N> &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const matrix<X,M,N> &A, const X &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const X &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const X &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const matrix<X,M,N> &A, const X &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const X &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const X &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const X &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const X &A, const matrix<X,M,N> &B)
{
  reg::matrix<X,M,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X, size_t M, size_t N>
inline std::vector<size_t> matrix<X,M,N>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline std::vector<size_t> matrix<X,M,N>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline size_t matrix<X,M,N>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X, size_t M, size_t N>
inline X matrix<X,M,N>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X, size_t M, size_t N>
inline X matrix<X,M,N>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X, size_t M, size_t N>
inline X matrix<X,M,N>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X, size_t M, size_t N>
inline double matrix<X,M,N>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X, size_t M, size_t N>
inline double matrix<X,M,N>::average(const matrix<X,M,N> &weights, bool norm) const
{
  assert( M == weights.shape(0) );
  assert( N == weights.shape(1) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}


// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t M, size_t N>
inline void matrix<X,M,N>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < M ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( j != N-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else            std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
inline std::ostream& operator<<(std::ostream& out, const matrix<X,M,N>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < M ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != N-1 ) out << ", ";
      else if ( i != M-1 ) out << ";" << std::endl;
      else                 out << ";";
    }
  }

  return out;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

