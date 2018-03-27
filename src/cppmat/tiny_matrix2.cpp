/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX_CPP
#define CPPMAT_TINY_MATRIX_CPP

#include "tiny_matrix2.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2()
{
  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2(X D)
{
  // copy input
  for ( size_t i = 0; i < m_size ; ++i ) m_container[i] = D;

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<typename Iterator>
inline matrix2<X,m,n>::matrix2(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // initialize counter
  size_t i = 0;

  // copy input
  for (auto it = first; it != last; ++it)
  {
    m_container[i] = *it; ++i;
  }

  // point to local data container
  m_data = &m_container[0];
}

// =================================================================================================
// copy constructor
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2(const matrix2<X,m,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// =================================================================================================
// assignment operator
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator= (const matrix2<X,m,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];

  // return pointer to current instance
  return *this;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::map(X *D)
{
  m_data = D;
}

// =================================================================================================
// copy from external pointer
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::copy(const X *D)
{
  std::copy(D, D+m_size, &m_container[0]);

  m_data = &m_container[0];
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t m, size_t n> inline size_t matrix2<X,m,n>::size() const { return m_size; }
template<class X, size_t m, size_t n> inline size_t matrix2<X,m,n>::ndim() const { return 2;      }

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

template<class X, size_t m, size_t n> inline X&       matrix2<X,m,n>::operator[](size_t i)       { return m_data[i]; }
template<class X, size_t m, size_t n> inline const X& matrix2<X,m,n>::operator[](size_t i) const { return m_data[i]; }

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline X& matrix2<X,m,n>::operator()(size_t a)
{
  return m_data[a*n];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix2<X,m,n>::operator()(size_t a) const
{
  return m_data[a*n];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline X& matrix2<X,m,n>::operator()(size_t a, size_t b)
{
  return m_data[a*n+b];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix2<X,m,n>::operator()(size_t a, size_t b) const
{
  return m_data[a*n+b];
}


// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X, size_t m, size_t n> inline X*       matrix2<X,m,n>::data()        { return &m_data[0];          }
template<class X, size_t m, size_t n> inline const X* matrix2<X,m,n>::data() const  { return &m_data[0];          }
template<class X, size_t m, size_t n> inline auto     matrix2<X,m,n>::begin()       { return &m_data[0];          }
template<class X, size_t m, size_t n> inline auto     matrix2<X,m,n>::begin() const { return &m_data[0];          }
template<class X, size_t m, size_t n> inline auto     matrix2<X,m,n>::end()         { return &m_data[0] + m_size; }
template<class X, size_t m, size_t n> inline auto     matrix2<X,m,n>::end() const   { return &m_data[0] + m_size; }

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::setConstant(X D)
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::setZero()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::setOnes()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::zeros()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::ones()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator*= (const matrix2<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator/= (const matrix2<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator+= (const matrix2<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator-= (const matrix2<X,m,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

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
inline std::ostream& operator<<(std::ostream& out, matrix2<X,m,n>& src)
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

#endif

