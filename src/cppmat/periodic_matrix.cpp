/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX2_CPP
#define CPPMAT_PERIODIC_MATRIX2_CPP

// -------------------------------------------------------------------------------------------------

#include "periodic_matrix2.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline matrix2<X>::matrix2(size_t m, size_t n)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(m,n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> matrix2<X>::Arange(size_t m, size_t n)
{
  // allocate matrix
  matrix2<X> out(m,n);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> matrix2<X>::Zero(size_t m, size_t n)
{
  // allocate matrix
  matrix2<X> out(m,n);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> matrix2<X>::Ones(size_t m, size_t n)
{
  // allocate matrix
  matrix2<X> out(m,n);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> matrix2<X>::Constant(size_t m, size_t n, X D)
{
  // allocate matrix
  matrix2<X> out(m,n);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline matrix2<X> matrix2<X>::Copy(size_t m, size_t n, Iterator first, Iterator last)
{
  // allocate matrix
  matrix2<X> out(m,n);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void matrix2<X>::resize(size_t m, size_t n)
{
  m_m    = m;
  m_n    = n;
  m_m_i  = static_cast<int>(m);
  m_n_i  = static_cast<int>(n);
  m_size = m*n;

  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::reshape(size_t m, size_t n)
{
  assert( m_size == m*n );

  resize(m,n);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t matrix2<X>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix2<X>::ndim() const
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix2<X>::shape(int i) const
{
  i = ( i < 0 ) ? i + 2 : ( i >= 2 ) ? i - 2 : i ;

  if ( i == 0 ) return m_m;
  if ( i == 1 ) return m_n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix2<X>::shape(size_t i) const
{
  if ( i == 0 ) return m_m;
  if ( i == 1 ) return m_n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix2<X>::shape() const
{
  std::vector<size_t> ret(2);

  ret[0] = m_m;
  ret[1] = m_n;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix2<X>::strides(bool bytes) const
{
  std::vector<size_t> ret(2);

  ret[0] = m_n;
  ret[1] = 1;

  if ( bytes )
    for ( size_t i = 0 ; i < 2 ; ++i )
      ret[i] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators
// =================================================================================================

template<class X>
inline X& matrix2<X>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix2<X>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------


template<class X>
inline X& matrix2<X>::operator()(int a)
{
  a = ( a < 0 ) ? a + m_m_i : ( a >= m_m_i ) ? a - m_m_i : a ;

  return m_data[a*m_n];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix2<X>::operator()(int a) const
{
  a = ( a < 0 ) ? a + m_m_i : ( a >= m_m_i ) ? a - m_m_i : a ;

  return m_data[a*m_n];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix2<X>::operator()(int a, int b)
{
  a = ( a < 0 ) ? a + m_m_i : ( a >= m_m_i ) ? a - m_m_i : a ;
  b = ( b < 0 ) ? b + m_n_i : ( b >= m_n_i ) ? b - m_n_i : b ;

  return m_data[a*m_n+b];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix2<X>::operator()(int a, int b) const
{
  a = ( a < 0 ) ? a + m_m_i : ( a >= m_m_i ) ? a - m_m_i : a ;
  b = ( b < 0 ) ? b + m_n_i : ( b >= m_n_i ) ? b - m_n_i : b ;

  return m_data[a*m_n+b];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline X* matrix2<X>::data()
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* matrix2<X>::data() const
{
  return m_data.data();
}

// =================================================================================================
// iterators
// =================================================================================================

template<class X>
inline auto matrix2<X>::begin()
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::end()
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::end() const
{
  return m_data.end();
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void matrix2<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void matrix2<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline matrix2<X>& matrix2<X>::operator*= (const matrix2<X> &B)
{
  assert( m_m == B.shape(0) );
  assert( m_n == B.shape(1) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator/= (const matrix2<X> &B)
{
  assert( m_m == B.shape(0) );
  assert( m_n == B.shape(1) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator+= (const matrix2<X> &B)
{
  assert( m_m == B.shape(0) );
  assert( m_n == B.shape(1) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator-= (const matrix2<X> &B)
{
  assert( m_m == B.shape(0) );
  assert( m_n == B.shape(1) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator* (const matrix2<X> &A, const matrix2<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator/ (const matrix2<X> &A, const matrix2<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator+ (const matrix2<X> &A, const matrix2<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator- (const matrix2<X> &A, const matrix2<X> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator* (const matrix2<X> &A, const X &B)
{
  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator/ (const matrix2<X> &A, const X &B)
{
  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator+ (const matrix2<X> &A, const X &B)
{
  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator- (const matrix2<X> &A, const X &B)
{
  matrix2<X> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator* (const X &A, const matrix2<X> &B)
{
  matrix2<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator/ (const X &A, const matrix2<X> &B)
{
  matrix2<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator+ (const X &A, const matrix2<X> &B)
{
  matrix2<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X> operator- (const X &A, const matrix2<X> &B)
{
  matrix2<X> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X matrix2<X>::minCoeff() const
{
  return *std::min_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X matrix2<X>::maxCoeff() const
{
  return *std::max_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X matrix2<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double matrix2<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double matrix2<X>::average(const matrix2<X> &weights) const
{
  assert( m_m == weights.shape(0) );
  assert( m_n == weights.shape(1) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void matrix2<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m_m ; ++i ) {
    for ( size_t j = 0 ; j < m_n ; ++j ) {
      if ( j != m_n-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else              std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const matrix2<X>& src)
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

#endif

