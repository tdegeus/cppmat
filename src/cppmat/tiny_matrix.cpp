/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX_CPP
#define CPPMAT_TINY_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "tiny_matrix.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix<X,m,n> matrix<X,m,n>::Arange()
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> matrix<X,m,n>::Zero()
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> matrix<X,m,n>::Ones()
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> matrix<X,m,n>::Constant(X D)
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<typename Iterator>
inline matrix<X,m,n> matrix<X,m,n>::Copy(Iterator first)
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setCopy(first,first+out.size());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<typename Iterator>
inline matrix<X,m,n> matrix<X,m,n>::Copy(Iterator first, Iterator last)
{
  // call basic constructor
  matrix<X,m,n> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::Size()
{
  return m_size;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::ndim() const
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::rows() const
{
  return m_m;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::cols() const
{
  return m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::shape(int i) const
{
  // check axis: (0,1) or (-1,-2)
  assert( i  <  2 );
  assert( i >= -2 );

  // correct periodic index
  i = ( 2 + (i%2) ) % 2;

  // return shape
  if ( i == 0 ) return m_m;
  else          return m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::shape(size_t i) const
{
  // check axis: (0,1)
  assert( i < 2 );

  // return shape
  if ( i == 0 ) return m_m;
  else          return m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix<X,m,n>::shape() const
{
  std::vector<size_t> ret(2);

  ret[0] = m_m;
  ret[1] = m_n;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix<X,m,n>::strides(bool bytes) const
{
  std::vector<size_t> ret(2);

  ret[0] = m_n;
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

template<class X, size_t m, size_t n>
inline X& matrix<X,m,n>::operator[](size_t i)
{
  assert( i < m_size );

  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix<X,m,n>::operator[](size_t i) const
{
  assert( i < m_size );

  return m_data[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X, size_t m, size_t n>
inline X& matrix<X,m,n>::operator()(size_t a)
{
  assert( a < m_m );

  return m_data[a*m_n];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix<X,m,n>::operator()(size_t a) const
{
  assert( a < m_m );

  return m_data[a*m_n];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline X& matrix<X,m,n>::operator()(size_t a, size_t b)
{
  assert( a < m_m );
  assert( b < m_n );

  return m_data[a*m_n+b];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix<X,m,n>::operator()(size_t a, size_t b) const
{
  assert( a < m_m );
  assert( b < m_n );

  return m_data[a*m_n+b];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X, size_t m, size_t n>
template<class Iterator>
inline X& matrix<X,m,n>::at(Iterator first, Iterator last)
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

  return m_data[a*m_n+b];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<class Iterator>
inline const X& matrix<X,m,n>::at(Iterator first, Iterator last) const
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

  return m_data[a*m_n+b];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::compress(size_t a) const
{
  assert( a < m_m );

  return a*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix<X,m,n>::compress(size_t a, size_t b) const
{
  assert( a < m_m );
  assert( b < m_n );

  return a*m_n+b;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix<X,m,n>::decompress(size_t i) const
{
  // check input
  assert( i < m_size );

  // allocate array-index
  std::vector<size_t> idx(2);

  // reconstruct
  idx[1] = i % m_n;
  idx[0] = ( i - idx[1] ) / m_n;

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t m, size_t n>
inline X* matrix<X,m,n>::data()
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X* matrix<X,m,n>::data() const
{
  return std::begin(m_data);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::begin()
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::begin() const
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::end()
{
  return std::begin(m_data) + m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::end() const
{
  return std::begin(m_data) + m_size;
}

// =================================================================================================
// iterators : beginRow() and endRow()
// =================================================================================================

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::beginRow(size_t a)
{
  assert( a < m_m );

  return std::begin(m_data) + a*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::beginRow(size_t a) const
{
  assert( a < m_m );

  return std::begin(m_data) + a*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::endRow(size_t a)
{
  assert( a < m_m );

  return std::begin(m_data) + (a+1)*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::endRow(size_t a) const
{
  assert( a < m_m );

  return std::begin(m_data) + (a+1)*m_n;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::index(size_t i)
{
  assert( i < m_size );

  return std::begin(m_data) + i;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::index(size_t i) const
{
  assert( i < m_size );

  return std::begin(m_data) + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::item(size_t a)
{
  assert( a < m_m );

  return std::begin(m_data) + a*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::item(size_t a) const
{
  assert( a < m_m );

  return std::begin(m_data) + a*m_n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::item(size_t a, size_t b)
{
  assert( a < m_m );
  assert( b < m_n );

  return std::begin(m_data) + a*m_n+b;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix<X,m,n>::item(size_t a, size_t b) const
{
  assert( a < m_m );
  assert( b < m_n );

  return std::begin(m_data) + a*m_n+b;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<class Iterator>
inline void matrix<X,m,n>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, std::begin(m_data));
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator*= (const matrix<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator/= (const matrix<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator+= (const matrix<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator-= (const matrix<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n>& matrix<X,m,n>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator* (const matrix<X,m,n> &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator/ (const matrix<X,m,n> &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator+ (const matrix<X,m,n> &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator- (const matrix<X,m,n> &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator* (const matrix<X,m,n> &A, const X &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator/ (const matrix<X,m,n> &A, const X &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator+ (const matrix<X,m,n> &A, const X &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator- (const matrix<X,m,n> &A, const X &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator* (const X &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator/ (const X &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator+ (const X &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix<X,m,n> operator- (const X &A, const matrix<X,m,n> &B)
{
  matrix<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix<X,m,n>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix<X,m,n>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X, size_t m, size_t n>
inline X matrix<X,m,n>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X, size_t m, size_t n>
inline X matrix<X,m,n>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X, size_t m, size_t n>
inline X matrix<X,m,n>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X, size_t m, size_t n>
inline double matrix<X,m,n>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X, size_t m, size_t n>
inline double matrix<X,m,n>::average(const matrix<X,m,n> &weights, bool norm) const
{
  assert( m_m == weights.shape(0) );
  assert( m_n == weights.shape(1) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// =================================================================================================
// basic algebra : absolute value
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::abs()
{
  for ( auto &i : m_data )
    i = std::abs(i);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix<X,m,n>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m_m ; ++i ) {
    for ( size_t j = 0 ; j < m_n ; ++j ) {
      if ( j != m_n-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else              std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::ostream& operator<<(std::ostream& out, const matrix<X,m,n>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < m ; ++i ) {
    for ( size_t j = 0 ; j < n ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != n-1 ) out << ", ";
      else if ( i != m-1 ) out << ";" << std::endl;
      else                 out << ";";
    }
  }

  return out;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

