/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VECTOR_H
#define CPPMAT_VECTOR_H

#include "macros.h"

namespace cppmat {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X>
class vector
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_n=0;     // number of columns
  size_t         m_size=0;  // total size

public:

  // constructors
  vector(){};
  vector(size_t n);
  vector(size_t n, X D);
  vector(size_t n, const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // resize
  void resize(size_t n);

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;

  // arithmetic operators
  vector<X>& operator*= (const vector<X> &B);
  vector<X>& operator/= (const vector<X> &B);
  vector<X>& operator+= (const vector<X> &B);
  vector<X>& operator-= (const vector<X> &B);
  vector<X>& operator*= (             X   B);
  vector<X>& operator/= (             X   B);
  vector<X>& operator+= (             X   B);
  vector<X>& operator-= (             X   B);

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
  double average(const vector<X> &weights) const;

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
  operator vector<U> () const;

  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator std::vector<U> () const;

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,1,Eigen::Dynamic,Eigen::RowMajor> () const;
  #endif

  #ifdef CPPMAT_EIGEN
  template<typename U,typename V=X,typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,Eigen::Dynamic,1,Eigen::ColMajor> () const;
  #endif

}; // class vector

template<class X> inline vector<X> operator* (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator/ (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator+ (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator- (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator* (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator/ (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator+ (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator- (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator* (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator/ (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator+ (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator- (const        X  &A, const vector<X> &B);

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline vector<X>::vector(size_t n)
{
  resize(n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector(size_t n, X D)
{
  resize(n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector(size_t n, const X *D)
{
  resize(n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D[i];
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t vector<X>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::ndim() const
{
  return 1;
}

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
// resize
// =================================================================================================

template<class X>
inline void vector<X>::resize(size_t n)
{
  m_n    = n;
  m_size = n;

  m_data.resize(m_size);
}

// =================================================================================================
// index operators
// =================================================================================================

template<class X>
inline X& vector<X>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& vector<X>::operator()(size_t a)
{
  return m_data[a];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator()(size_t a) const
{
  return m_data[a];
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
};

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
};

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  assert( m_n == B.shape(0) );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
};

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
inline vector<X>& vector<X>::operator*= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (X B)
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
// pointers / iterators
// =================================================================================================

template<class X>
inline X* vector<X>::data()
{
  return m_data.data ();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::begin()
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::end()
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* vector<X>::data() const
{
  return m_data.data ();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::end() const
{
  return m_data.end();
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
// basic initialization
// =================================================================================================

template<class X>
inline void vector<X>::setConstant(X D)
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setZero()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setOnes()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::zeros()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::ones()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// conversion operators
// =================================================================================================

template<class X>
template<typename U, typename V, typename E>
inline vector<X>::operator vector<U> () const
{
  vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = static_cast<U>( m_data[i] );

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename U, typename V, typename E>
inline vector<X>::operator std::vector<U> () const
{
  std::vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = static_cast<U>( m_data[i] );

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
inline vector<X>::operator Eigen::Matrix<U,1,Eigen::Dynamic,Eigen::RowMajor> () const
{
  Eigen::Matrix<U,1,Eigen::Dynamic,Eigen::RowMajor> out(m_n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X>
template<typename U, typename V, typename E>
inline vector<X>::operator Eigen::Matrix<U,Eigen::Dynamic,1,Eigen::ColMajor> () const
{
  Eigen::Matrix<U,Eigen::Dynamic,1,Eigen::ColMajor> out(m_n);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  std::vector<size_t> s = strides();

  for ( size_t h = 0 ; h < shape(0)-1 ; ++h )
    std::printf((fmt+",").c_str(),m_data[h]);

  std::printf((fmt+"\n").c_str(),m_data[shape(0)-1]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

// =================================================================================================

} // namespace ...

#endif

