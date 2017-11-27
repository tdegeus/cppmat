/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_VECTOR_H
#define CPPMAT_TINY_VECTOR_H

#include "macros.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::vector
// =================================================================================================

template<class X, size_t n>
class vector
{
private:

  X      m_container[n];    // data container
  X     *m_data;            // pointer to container (may point outside)
  size_t m_size=n;          // total number of entries

public:

  // constructors
  vector();
  vector(X D);
  vector(const X *D);
  vector(const vector<X,n> &D);

  // copy constructor
  vector<X,n>& operator= (const vector<X,n> &D);

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

  // arithmetic operators
  vector<X,n>& operator*= (const vector<X,n> &B);
  vector<X,n>& operator/= (const vector<X,n> &B);
  vector<X,n>& operator+= (const vector<X,n> &B);
  vector<X,n>& operator-= (const vector<X,n> &B);
  vector<X,n>& operator*= (             X     B);
  vector<X,n>& operator/= (             X     B);
  vector<X,n>& operator+= (             X     B);
  vector<X,n>& operator-= (             X     B);

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
  double average(const vector<X,n> &weights) const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // formatted print; NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // conversion operators
  template<typename U, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator vector<U,n> () const;

  template<typename U, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator std::vector<U> () const;

  #ifdef CPPMAT_EIGEN
  template<typename U, int N, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,1,N,Eigen::RowMajor> () const;
  #endif

  #ifdef CPPMAT_EIGEN
  template<typename U, int N, typename V=X, typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator Eigen::Matrix<U,N,1,Eigen::ColMajor> () const;
  #endif

}; // class vector

// arithmetic operators
template<class X, size_t n> inline vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator* (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator/ (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator+ (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator- (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator* (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator/ (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator+ (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator- (const        X    &A, const vector<X,n> &B);


// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>::vector()
{
  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>::vector(X D)
{
  // copy input
  for ( size_t i = 0; i < m_size ; ++i ) m_container[i] = D;

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>::vector(const X *D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>::vector(const vector<X,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// =================================================================================================
// copy constructors
// =================================================================================================

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator= (const vector<X,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_container[i] = D[i];
  // point to local data container
  m_data = &m_container[0];
  // return pointer to current instance
  return *this;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::map(X *D)
{
  m_data = D;
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
inline vector<X,n>& vector<X,n>::operator*= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator/= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator+= (X B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline vector<X,n>& vector<X,n>::operator-= (X B)
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
// pointers / iterators
// =================================================================================================

template<class X, size_t n>
inline X* vector<X,n>::data()
{
  return &m_data[0];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline auto vector<X,n>::begin()
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
inline const X* vector<X,n>::data() const
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
inline auto vector<X,n>::end() const
{
  return &m_data[0] + m_size;
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
// basic initialization
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::setConstant(X D)
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setZero()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::setOnes()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::zeros()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline void vector<X,n>::ones()
{
  for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// conversion operators
// =================================================================================================

template<class X, size_t n>
template<typename U, typename V, typename E>
inline vector<X,n>::operator vector<U,n> () const
{
  vector<U,n> out;

  for ( size_t i = 0 ; i < n ; ++i )
    out[i] = static_cast<U>(m_data[i]);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
template<typename U, typename V, typename E>
inline vector<X,n>::operator std::vector<U> () const
{
  std::vector<U> out(m_size);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out[i] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X, size_t n>
template<typename U, int N, typename V, typename E>
inline vector<X,n>::operator Eigen::Matrix<U,1,N,Eigen::RowMajor> () const
{
  Eigen::Matrix<U,1,N,Eigen::RowMajor> out;

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifdef CPPMAT_EIGEN
template<class X, size_t n>
template<typename U, int N, typename V, typename E>
inline vector<X,n>::operator Eigen::Matrix<U,N,1,Eigen::ColMajor> () const
{
  Eigen::Matrix<U,N,1,Eigen::ColMajor> out;

  for ( size_t i = 0 ; i < m_size ; ++i )
    out(i) = m_data[i];

  return out;
}
#endif

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t n>
inline void vector<X,n>::printf(std::string fmt) const
{
  std::vector<size_t> s = strides();

  for ( size_t h = 0 ; h < shape(0)-1 ; ++h )
    std::printf((fmt+",").c_str(),m_data[h]);

  std::printf((fmt+"\n").c_str(),m_data[shape(0)-1]);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t n>
inline std::ostream& operator<<(std::ostream& out, vector<X,n>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

// =================================================================================================

}} // namespace ...

#endif

