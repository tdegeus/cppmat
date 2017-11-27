/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX2_H
#define CPPMAT_TINY_MATRIX2_H

#include "macros.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::matrix2
// =================================================================================================

template<class X, size_t m, size_t n>
class matrix2
{
private:

  X      m_container[m*n];  // data container
  X     *m_data;            // pointer to container (may point outside)
  size_t m_size=m*n;        // total number of entries

public:

  // constructors
  matrix2();
  matrix2(X D);
  matrix2(const X *D);
  matrix2(const matrix2<X,m,n> &D);

  // copy constructor
  matrix2<X,m,n>& operator= (const matrix2<X,m,n> &D);

  // map external pointer
  void map(X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;

  // arithmetic operators
  matrix2<X,m,n>& operator*= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator/= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator+= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator-= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator*= (              X       B);
  matrix2<X,m,n>& operator/= (              X       B);
  matrix2<X,m,n>& operator+= (              X       B);
  matrix2<X,m,n>& operator-= (              X       B);

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const matrix2<X,m,n> &weights) const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // formatted print; NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // conversion operators
  template<typename U, size_t M, size_t N, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator matrix2<U,M,N> () const;

  template<typename U, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator std::vector<U> () const;

  #ifdef CPPMAT_EIGEN
  template<typename U, int M, int N, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,M,N,Eigen::RowMajor> () const;
  #endif

  #ifdef CPPMAT_EIGEN
  template<typename U, int M, int N, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,M,N,Eigen::ColMajor> () const;
  #endif

}; // class matrix2

// arithmetic operators
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator* (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator/ (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator+ (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline matrix2<X,m,n> operator- (const         X      &A, const matrix2<X,m,n> &B);

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
inline matrix2<X,m,n>::matrix2(const X *D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>::matrix2(const matrix2<X,m,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// =================================================================================================
// copy constructors
// =================================================================================================

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator= (const matrix2<X,m,n> &D)
{
  // - copy input
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_container[i] = D[i];
  // - point to local data container
  m_data = &m_container[0];
  // - return pointer to current instance
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
// get dimensions
// =================================================================================================

template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::size() const
{
  return m_size;

// -------------------------------------------------------------------------------------------------
}
template<class X, size_t m, size_t n>
inline size_t matrix2<X,m,n>::ndim() const
{
  return 2;
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

template<class X, size_t m, size_t n>
inline X& matrix2<X,m,n>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X& matrix2<X,m,n>::operator[](size_t i) const
{
  return m_data[i];
}

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
inline matrix2<X,m,n>& matrix2<X,m,n>::operator*= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator/= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator+= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline matrix2<X,m,n>& matrix2<X,m,n>::operator-= (X B)
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
// pointers / iterators
// =================================================================================================

template<class X, size_t m, size_t n>
inline X* matrix2<X,m,n>::data()
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix2<X,m,n>::begin()
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix2<X,m,n>::end()
{
  return &m_data[0] + m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline const X* matrix2<X,m,n>::data() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix2<X,m,n>::begin() const
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline auto matrix2<X,m,n>::end() const
{
  return &m_data[0] + m_size;
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
// conversion operators
// =================================================================================================

template<class X, size_t m, size_t n>
template<typename U, size_t M, size_t N, typename V, typename E>
inline matrix2<X,m,n>::operator matrix2<U,M,N> () const
{
  matrix2<U,M,N> out;

  for ( size_t i = 0 ; i < M*N ; ++i )
    out[i] = static_cast<U>(m_data[i]);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
template<typename U, typename V, typename E>
inline matrix2<X,m,n>::operator std::vector<U> () const
{
  std::vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X, size_t m, size_t n>
template<typename U, int M, int N, typename V, typename E>
inline matrix2<X,m,n>::operator Eigen::Matrix<U,M,N,Eigen::RowMajor> () const
{
  Eigen::Matrix<U,M,N,Eigen::RowMajor> out;

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X, size_t m, size_t n>
template<typename U, int M, int N, typename V, typename E>
inline matrix2<X,m,n>::operator Eigen::Matrix<U,M,N,Eigen::ColMajor> () const
{
  Eigen::Matrix<U,M,N,Eigen::RowMajor> out;

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t m, size_t n>
inline void matrix2<X,m,n>::printf(std::string fmt) const
{
  std::vector<size_t> s = strides();

  for ( size_t h = 0 ; h < shape(0) ; ++h )
  {
    for ( size_t i = 0 ; i < shape(1)-1 ; ++i )
      std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]]);

    std::printf((fmt+";\n").c_str(),m_data[h*s[0]+(shape(1)-1)*s[1]]);
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t m, size_t n>
inline std::ostream& operator<<(std::ostream& out, matrix2<X,m,n>& src)
{
  for ( size_t i = 0 ; i < src.shape(0) ; ++i )
  {
    for ( size_t j = 0 ; j < src.shape(1)-1 ; ++j )
      out << src(i,j) << ", ";

    out << src(i,src.shape(1)-1) << "; " << std::endl;
  }

  return out;
}

// =================================================================================================

}} // namespace ...

#endif

