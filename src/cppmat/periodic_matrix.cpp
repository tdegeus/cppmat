/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX_CPP
#define CPPMAT_PERIODIC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "periodic_matrix.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline matrix<X>::matrix(size_t m, size_t n)
{
  // store shape, and other size parameters, allocate "mData"
  resize(m,n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> matrix<X>::Arange(size_t m, size_t n)
{
  // call basic constructor
  matrix<X> out(m,n);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> matrix<X>::Zero(size_t m, size_t n)
{
  // call basic constructor
  matrix<X> out(m,n);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> matrix<X>::Ones(size_t m, size_t n)
{
  // call basic constructor
  matrix<X> out(m,n);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> matrix<X>::Constant(size_t m, size_t n, X D)
{
  // call basic constructor
  matrix<X> out(m,n);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline matrix<X> matrix<X>::Copy(size_t m, size_t n, Iterator first, Iterator last)
{
  // call basic constructor
  matrix<X> out(m,n);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void matrix<X>::resize(size_t m, size_t n)
{
  // update shape/size
  M     = m;
  N     = n;
  MI    = static_cast<int>(m);
  NI    = static_cast<int>(n);
  mSize = m*n;

  // allocate data
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::reshape(size_t m, size_t n)
{
  // check that the size is unchanged
  assert( mSize == m*n );

  // process new shape
  resize(m,n);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t matrix<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::ndim() const
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::rows() const
{
  return M;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::cols() const
{
  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::shape(int i) const
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

template<class X>
inline size_t matrix<X>::shape(size_t i) const
{
  // check axis: (0,1)
  assert( i < 2 );

  // return shape
  if ( i == 0 ) return M;
  else          return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix<X>::shape() const
{
  std::vector<size_t> ret(2);

  ret[0] = M;
  ret[1] = N;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix<X>::strides(bool bytes) const
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

template<class X>
inline X& matrix<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline X& matrix<X>::operator()(int a)
{
  a = ( MI + (a % MI) ) % MI;

  assert( a < MI );

  return mData[a*N];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a) const
{
  a = ( MI + (a % MI) ) % MI;

  assert( a < MI );

  return mData[a*N];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b)
{
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  assert( a < MI );
  assert( b < NI );

  return mData[a*N+b];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b) const
{
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  assert( a < MI );
  assert( b < NI );

  return mData[a*N+b];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X>
template<class Iterator>
inline X& matrix<X>::at(Iterator first, Iterator last)
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

  // correct for periodicity
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  return mData[a*N+b];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline const X& matrix<X>::at(Iterator first, Iterator last) const
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

  // correct for periodicity
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  return mData[a*N+b];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline size_t matrix<X>::compress(int a) const
{
  a = ( MI + (a % MI) ) % MI;

  assert( a < M );

  return static_cast<size_t>(a*NI);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::compress(int a, int b) const
{
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  assert( a < M );
  assert( b < N );

  return static_cast<size_t>(a*NI+b);
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X>
inline std::vector<size_t> matrix<X>::decompress(size_t i) const
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

template<class X>
inline X* matrix<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* matrix<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline auto matrix<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : beginRow() and endRow()
// =================================================================================================

template<class X>
inline auto matrix<X>::beginRow(size_t a)
{
  assert( a < M );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::beginRow(size_t a) const
{
  assert( a < M );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::endRow(size_t a)
{
  assert( a < M );

  return begin() + (a+1)*N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::endRow(size_t a) const
{
  assert( a < M );

  return begin() + (a+1)*N;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto matrix<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline auto matrix<X>::item(int a)
{
  a = ( MI + (a % MI) ) % MI;

  assert( a < MI );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::item(int a) const
{
  a = ( MI + (a % MI) ) % MI;

  assert( a < MI );

  return begin() + a*N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::item(int a, int b)
{
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  assert( a < MI );
  assert( b < NI );

  return begin() + a*N+b;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::item(int a, int b) const
{
  a = ( MI + (a % MI) ) % MI;
  b = ( NI + (b % NI) ) % NI;

  assert( a < MI );
  assert( b < NI );

  return begin() + a*N+b;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void matrix<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void matrix<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last - first );

  std::copy(first, last, begin());
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline matrix<X>& matrix<X>::operator*= (const matrix<X> &B)
{
  assert( M == B.shape(0) );
  assert( N == B.shape(1) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator/= (const matrix<X> &B)
{
  assert( M == B.shape(0) );
  assert( N == B.shape(1) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator+= (const matrix<X> &B)
{
  assert( M == B.shape(0) );
  assert( N == B.shape(1) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator-= (const matrix<X> &B)
{
  assert( M == B.shape(0) );
  assert( N == B.shape(1) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X>
inline std::vector<size_t> matrix<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X>
inline X matrix<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::minCoeff(size_t axis) const
{
  // check input
  assert( axis < 2 );

  // contraction of columns
  if ( axis == 1 )
  {
    vector<X> out = vector<X>::Constant(M, this->maxCoeff());

    for ( size_t i = 0 ; i < M ; ++i )
      for ( size_t j = 0 ; j < N ; ++j )
        out[i] = std::min( out[i], mData[i*N+j] );

    return out;
  }

  // contraction of rows
  vector<X> out = vector<X>::Constant(N, this->maxCoeff());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      out[j] = std::min( out[j], mData[i*N+j] );

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::minCoeff(int axis) const
{
  // check axis: (0,1) or (-1,-2)
  assert( axis  <  2 );
  assert( axis >= -2 );

  // correct periodic axis
  axis = ( 2 + (axis%2) ) % 2;

  // compute
  return this->minCoeff(static_cast<size_t>(axis));
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X>
inline X matrix<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::maxCoeff(size_t axis) const
{
  // check input
  assert( axis < 2 );

  // contraction of columns
  if ( axis == 1 )
  {
    vector<X> out = vector<X>::Constant(M, this->minCoeff());

    for ( size_t i = 0 ; i < M ; ++i )
      for ( size_t j = 0 ; j < N ; ++j )
        out[i] = std::max( out[i], mData[i*N+j] );

    return out;
  }

  // contraction of rows
  vector<X> out = vector<X>::Constant(N, this->minCoeff());

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      out[j] = std::max( out[j], mData[i*N+j] );

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::maxCoeff(int axis) const
{
  // check axis: (0,1) or (-1,-2)
  assert( axis  <  2 );
  assert( axis >= -2 );

  // correct periodic axis
  axis = ( 2 + (axis%2) ) % 2;

  // compute
  return this->maxCoeff(static_cast<size_t>(axis));
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X>
inline X matrix<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::sum(size_t axis) const
{
  // check input
  assert( axis < 2 );

  // contraction of columns
  if ( axis == 1 )
  {
    vector<X> out = vector<X>::Zero(M);

    for ( size_t i = 0 ; i < M ; ++i )
      for ( size_t j = 0 ; j < N ; ++j )
        out[i] += mData[i*N+j];

    return out;
  }

  // contraction of rows
  vector<X> out = vector<X>::Zero(N);

  for ( size_t i = 0 ; i < M ; ++i )
    for ( size_t j = 0 ; j < N ; ++j )
      out[j] += mData[i*N+j];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::sum(int axis) const
{
  // check axis: (0,1) or (-1,-2)
  assert( axis  <  2 );
  assert( axis >= -2 );

  // correct periodic axis
  axis = ( 2 + (axis%2) ) % 2;

  // compute
  return this->sum(static_cast<size_t>(axis));
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X>
inline double matrix<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::mean(size_t axis) const
{
  matrix<X> weights = matrix<X>::Ones(M, N);

  return (*this).sum(axis) / weights.sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::mean(int axis) const
{
  matrix<X> weights = matrix<X>::Ones(M, N);

  return (*this).sum(axis) / weights.sum(axis);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X>
inline double matrix<X>::average(const matrix<X> &weights, bool norm) const
{
  assert( M == weights.shape(0) );
  assert( N == weights.shape(1) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::average(const matrix<X> &weights, size_t axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> matrix<X>::average(const matrix<X> &weights, int axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// =================================================================================================
// basic algebra : absolute value
// =================================================================================================

template<class X>
inline void matrix<X>::abs()
{
  for ( auto &i : mData )
    i = std::abs(i);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void matrix<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < M ; ++i ) {
    for ( size_t j = 0 ; j < N ; ++j ) {
      if ( j != N-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else            std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const matrix<X>& src)
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

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

