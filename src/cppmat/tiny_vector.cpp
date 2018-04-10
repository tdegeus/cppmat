/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_VECTOR_CPP
#define CPPMAT_TINY_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "tiny_vector.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>::vector(X D)
{
  // set all entries constant
  for ( size_t i = 0; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<typename Iterator>
inline vector<X,n>::vector(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // initialize counter
  size_t i = 0;

  // copy input
  for (auto it = first; it != last; ++it)
  {
    m_data[i] = *it; ++i;
  }

}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t n>
inline size_t vector<X,n>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::ndim() const
{
  return 1;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::shape(int i) const
{
  i = ( i < 0 ) ? i + 1 : ( i >= 1 ) ? i - 1 : i ;

  if ( i == 0 ) return n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::shape(size_t i) const
{
  if ( i == 0 ) return n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline std::vector<size_t> vector<X,n>::shape() const
{
  std::vector<size_t> ret(1);

  ret[0] = n;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline std::vector<size_t> vector<X,n>::strides(bool bytes) const
{
  std::vector<size_t> ret(1);

  ret[0] = 1;

  if ( bytes )
    ret[0] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators
// =================================================================================================

template<class X, size_t n>
inline X& vector<X,n>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X& vector<X,n>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline X& vector<X,n>::operator()(size_t a)
{
  return m_data[a];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X& vector<X,n>::operator()(size_t a) const
{
  return m_data[a];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t n>
inline X* vector<X,n>::data()
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X* vector<X,n>::data() const
{
  return &m_data[0];
}

// =================================================================================================
// iterators
// =================================================================================================

template<class X, size_t n>
inline auto vector<X,n>::begin()
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::begin() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::end()
{
  return &m_data[0] + m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::end() const
{
  return &m_data[0] + m_size;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::zeros()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::ones()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator*= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator/= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator+= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator-= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator* (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator/ (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator+ (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> operator- (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X, size_t n>
inline X vector<X,n>::min() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline X vector<X,n>::max() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline X vector<X,n>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline double vector<X,n>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m>
inline double vector<X,m>::average(const vector<X,m> &weights) const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < n ; ++j ) {
    if ( j != n-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else            std::printf((fmt + ";\n").c_str(), (*this)(j));
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline std::ostream& operator<<(std::ostream& out, vector<X,n>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < n ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src(j);
    if ( j != n-1 ) out << ", ";
  }

  return out;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

