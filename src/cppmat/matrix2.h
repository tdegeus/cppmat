/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MATRIX2_H
#define CPPMAT_MATRIX2_H

#include "macros.h"

namespace cppmat {

// =================================================================================================
// cppmat::matrix2
// =================================================================================================

template<class X>
class matrix2
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_m=0;     // number of rows
  size_t         m_n=0;     // number of columns
  size_t         m_size=0;  // total size

public:

  // constructors
  matrix2(){};
  matrix2(size_t m, size_t n);
  matrix2(size_t m, size_t n, X D);
  matrix2(size_t m, size_t n, const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;

  // arithmetic operators
  matrix2<X>& operator*= (const matrix2<X> &B);
  matrix2<X>& operator/= (const matrix2<X> &B);
  matrix2<X>& operator+= (const matrix2<X> &B);
  matrix2<X>& operator-= (const matrix2<X> &B);
  matrix2<X>& operator*= (X B);
  matrix2<X>& operator/= (X B);
  matrix2<X>& operator+= (X B);
  matrix2<X>& operator-= (X B);

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
  double average(const matrix2<X> &weights) const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // formatted print; NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // conversion operators
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator matrix2<U> () const;

  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator std::vector<U> () const;

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> () const;
  #endif

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> () const;
  #endif

}; // class matrix2

template<class X> inline matrix2<X> operator* (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator/ (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator+ (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator- (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator* (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator/ (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator+ (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator- (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator* (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator/ (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator+ (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator- (const         X  &A, const matrix2<X> &B);

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline matrix2<X>::matrix2(size_t m, size_t n)
{
  resize(m,n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>::matrix2(size_t m, size_t n, X D)
{
  resize(m,n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>::matrix2(size_t m, size_t n, const X *D)
{
  resize(m,n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D[i];
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
// resize
// =================================================================================================

template<class X>
inline void matrix2<X>::resize(size_t m, size_t n)
{
  m_m    = m;
  m_n    = n;
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
inline X& matrix2<X>::operator()(size_t a, size_t b)
{
  return m_data[a*m_n+b];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix2<X>::operator()(size_t a, size_t b) const
{
  return m_data[a*m_n+b];
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
inline matrix2<X>& matrix2<X>::operator*= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator/= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator+= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix2<X>& matrix2<X>::operator-= (X B)
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
// pointers / iterators
// =================================================================================================

template<class X>
inline X* matrix2<X>::data()
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::begin()
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
inline const X* matrix2<X>::data() const
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix2<X>::end() const
{
  return m_data.end();
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X matrix2<X>::min() const
{
  return *std::min_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X matrix2<X>::max() const
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
// basic initialization
// =================================================================================================

template<class X>
inline void matrix2<X>::setConstant(X D)
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::setZero()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::setOnes()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::zeros()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix2<X>::ones()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// conversion operators
// =================================================================================================

template<class X>
template<typename U, typename V, typename E>
inline matrix2<X>::operator matrix2<U> () const
{
  matrix2<U> out(shape(0),shape(1));

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = static_cast<U>( m_data[i] );

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U, typename V, typename E>
inline matrix2<X>::operator std::vector<U> () const
{
  std::vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
inline matrix2<X>::operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> () const
{
  Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> out(m_m,m_n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
inline matrix2<X>::operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> () const
{
  Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> out(m_m,m_n);

  for ( size_t i = 0 ; i < m_m ; ++i )
    for ( size_t j = 0 ; j < m_n ; ++j )
      out(i,j) = m_data[i*m_n+j];

  return out;
}
#endif

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void matrix2<X>::printf(std::string fmt) const
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

template<class X>
inline std::ostream& operator<<(std::ostream& out, matrix2<X>& src)
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

} // namespace cppmat

#endif

