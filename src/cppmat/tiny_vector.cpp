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
inline vector<X,n> vector<X,n>::Arange()
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> vector<X,n>::Zero()
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> vector<X,n>::Ones()
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n> vector<X,n>::Constant(X D)
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<typename Iterator>
inline vector<X,n> vector<X,n>::Copy(Iterator first)
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<typename Iterator>
inline vector<X,n> vector<X,n>::Copy(Iterator first, Iterator last)
{
  // call basic constructor
  vector<X,n> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::Size()
{
  return n;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t n>
inline size_t vector<X,n>::size() const
{
  return n;
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
  // check axis: (0) or (-1)
  assert( i  <  1 );
  assert( i >= -1 );

  return n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::shape(size_t i) const
{
  // check axis: (0)
  assert( i < 1 );

  return n;

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
// index operators : operator[...]
// =================================================================================================

template<class X, size_t n>
inline X& vector<X,n>::operator[](size_t i)
{
  assert( i < n );

  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X& vector<X,n>::operator[](size_t i) const
{
  assert( i < n );

  return m_data[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X, size_t n>
inline X& vector<X,n>::operator()(size_t a)
{
  assert( a < n );

  return m_data[a];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X& vector<X,n>::operator()(size_t a) const
{
  assert( a < n );

  return m_data[a];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t n>
inline X* vector<X,n>::data()
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline const X* vector<X,n>::data() const
{
  return std::begin(m_data);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t n>
inline auto vector<X,n>::begin()
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::begin() const
{
  return std::begin(m_data);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::end()
{
  return std::begin(m_data) + n;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::end() const
{
  return std::begin(m_data) + n;
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X, size_t n>
inline auto vector<X,n>::index(size_t i)
{
  assert( i < n );

  return std::begin(m_data) + i;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::index(size_t i) const
{
  assert( i < n );

  return std::begin(m_data) + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t n>
inline auto vector<X,n>::item(size_t a)
{
  assert( a < n );

  return std::begin(m_data) + a;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::item(size_t a) const
{
  assert( a < n );

  return std::begin(m_data) + a;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::setArange()
{
  for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setZero()
{
  for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setOnes()
{
  for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setConstant(X D)
{
  for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<class Iterator>
inline void vector<X,n>::setCopy(Iterator first)
{
  // copy
  std::copy(first, first+n, std::begin(m_data));
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<class Iterator>
inline void vector<X,n>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( n == last - first );

  // copy
  std::copy(first, last, std::begin(m_data));
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator*= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator/= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator+= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator-= (const vector<X,n> &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < n ; ++i )
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
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X, size_t n>
inline size_t vector<X,n>::argmin() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline size_t vector<X,n>::argmax() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X, size_t n>
inline X vector<X,n>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X, size_t n>
inline X vector<X,n>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X, size_t n>
inline X vector<X,n>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < n ; ++i )
    out += m_data[i];

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X, size_t n>
inline double vector<X,n>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(n);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X, size_t n>
inline double vector<X,n>::average(const vector<X,n> &weights, bool norm) const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < n ; ++i )
    out += m_data[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// =================================================================================================
// find all non-zero entries
// =================================================================================================

template<class X, size_t n>
inline cppmat::vector<size_t> vector<X,n>::where() const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < n ; ++i )
    if ( m_data[i] )
      ++nnz;

  cppmat::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < n ; ++i ) {
    if ( m_data[i] ) {
      out[j] = i;
      ++j;
    }
  }

  return out;
}

// =================================================================================================
// basic algebra : absolute value
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::abs()
{
  for ( auto &i : m_data )
    i = std::abs(i);
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
inline std::ostream& operator<<(std::ostream& out, const vector<X,n>& src)
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

