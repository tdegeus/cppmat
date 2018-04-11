/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_MATRIX_CPP
#define CPPMAT_VIEW_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "view_matrix2.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2()
{
  m_data = nullptr;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2(const X *D)
{
  m_data = D;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::map(const X *D)
{
  m_data = D;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::ndim() const
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::shape(int i) const
{
  i = ( i < 0 ) ? i + 2 : ( i >= 2 ) ? i - 2 : i ;

  if ( i == 0 ) return m;
  if ( i == 1 ) return n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::shape(size_t i) const
{
  if ( i == 0 ) return m;
  if ( i == 1 ) return n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix2<X,m,n>::shape() const
{
  std::vector<size_t> ret(2);

  ret[0] = m;
  ret[1] = n;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::vector<size_t> matrix2<X,m,n>::strides(bool bytes) const
{
  std::vector<size_t> ret(2);

  ret[0] = n;
  ret[1] = 1;

  if ( bytes )
    for ( size_t i = 0 ; i < 2 ; ++i )
      ret[i] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators
// =================================================================================================

template<class X, size_t m, size_t n> inline const X& matrix2<X,m,n>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix2<X,m,n>::operator()(size_t a) const
{
  return m_data[a*n];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix2<X,m,n>::operator()(size_t a, size_t b) const
{
  return m_data[a*n+b];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t m, size_t n> inline const X* matrix2<X,m,n>::data() const
{
  return &m_data[0];
}

// =================================================================================================
// iterators
// =================================================================================================

template<class X, size_t m, size_t n> inline auto matrix2<X,m,n>::begin() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n> inline auto matrix2<X,m,n>::end() const
{
  return &m_data[0] + m_size;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const X &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const X &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const X &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const X &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator* (const X &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator/ (const X &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator+ (const X &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline reg::matrix2<X,m,n> operator- (const X &A, const matrix2<X,m,n> &B)
{
  reg::matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X, size_t m, size_t n>
inline X matrix2<X,m,n>::min() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline X matrix2<X,m,n>::max() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline X matrix2<X,m,n>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline double matrix2<X,m,n>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline double matrix2<X,m,n>::average(const matrix2<X,m,n> &weights) const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m ; ++i ) {
    for ( size_t j = 0 ; j < n ; ++j ) {
      if ( j != n-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else            std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::ostream& operator<<(std::ostream& out, const matrix2<X,m,n>& src)
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

