/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_VECTOR_CPP
#define CPPMAT_PERIODIC_VECTOR_CPP

#include "periodic_vector.h"

namespace cppmat {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline vector<X>::vector(size_t n)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector(size_t n, X D)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(n);

  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline vector<X>::vector(Iterator first, Iterator last)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(last - first);

  // initialize counter
  size_t i = 0;

  // copy input
  for (auto it = first; it != last; ++it)
  {
    m_data[i] = *it; ++i;
  }
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void vector<X>::resize(size_t n)
{
  m_n    = n;
  m_n_i  = static_cast<int>(n);
  m_size = n;

  m_data.resize(m_size);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X> inline size_t vector<X>::size() const { return m_size; }
template<class X> inline size_t vector<X>::ndim() const { return 1;      }

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::shape(size_t i) const
{
  if ( i == 0 ) return m_n;

  assert( false );
  return 0;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> ret(1);

  ret[0] = m_n;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::strides(bool bytes) const
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

template<class X> inline X&       vector<X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& vector<X>::operator[](size_t i) const { return m_data[i]; }

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& vector<X>::operator()(int a)
{
  a = ( a < 0 ) ? a + m_n_i : ( a >= m_n_i ) ? a - m_n_i : a ;

  return m_data[a];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator()(int a) const
{
  a = ( a < 0 ) ? a + m_n_i : ( a >= m_n_i ) ? a - m_n_i : a ;

  return m_data[a];
}

// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X> inline X*       vector<X>::data()        { return m_data.data();  }
template<class X> inline const X* vector<X>::data() const  { return m_data.data();  }
template<class X> inline auto     vector<X>::begin()       { return m_data.begin(); }
template<class X> inline auto     vector<X>::begin() const { return m_data.begin(); }
template<class X> inline auto     vector<X>::end()         { return m_data.end();   }
template<class X> inline auto     vector<X>::end() const   { return m_data.end();   }

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void vector<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::zeros()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::ones()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline vector<X>& vector<X>::operator*= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X vector<X>::min() const
{
  return *std::min_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::max() const
{
  return *std::max_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double vector<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double vector<X>::average(const vector<X> &weights) const
{
  assert( m_n == weights.shape(0) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < m_n ; ++j ) {
    if ( j != m_n-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else              std::printf((fmt + ";\n").c_str(), (*this)(j));
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < src.size() ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src(j);
    if ( j != src.size()-1 ) out << ", ";
  }

  return out;
}

// =================================================================================================

}} // namespace ...

#endif

