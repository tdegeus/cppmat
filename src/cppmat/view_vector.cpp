/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_VECTOR_CPP
#define CPPMAT_VIEW_VECTOR_CPP

#include "view_vector.h"

namespace cppmat {
namespace view {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>::vector()
{
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>::vector(const X *D)
{
  m_data = D;
}

// =================================================================================================
// assignment operator
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator= (const vector<X,n> &D)
{
  m_data = &D[0];

  return *this;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::map(const X *D)
{
  m_data = D;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t n> inline size_t vector<X,n>::size() const { return m_size; }
template<class X, size_t n> inline size_t vector<X,n>::ndim() const { return 1;      }

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

template<class X, size_t n> inline const X& vector<X,n>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X& vector<X,n>::operator()(size_t a) const
{
  return m_data[a];
}

// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X, size_t n> inline const X* vector<X,n>::data() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n> inline auto vector<X,n>::begin() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n> inline auto vector<X,n>::end() const
{
  return &m_data[0] + m_size;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator* (const vector<X,n> &A, const X &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator/ (const vector<X,n> &A, const X &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator+ (const vector<X,n> &A, const X &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator- (const vector<X,n> &A, const X &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator* (const X &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator/ (const X &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator+ (const X &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline cppmat::tiny::vector<X,n> operator- (const X &A, const vector<X,n> &B)
{
  cppmat::tiny::vector<X,n> C;

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

#endif

