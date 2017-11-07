/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX2_H
#define CPPMAT_PERIODIC_MATRIX2_H

#include "macros.h"

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::matrix2
// =================================================================================================

template <class T> class matrix2
{
private:

  std::vector<T> m_data;   // data array
  size_t         m_m=0;    // number of rows
  size_t         m_n=0;    // number of columns

public:

  // (copy) constructor
  // ------------------

  matrix2               (const matrix2<T> &) = default;
  matrix2<T>& operator= (const matrix2<T> &) = default;
  matrix2<T>(){};

  // explicit constructors
  // ---------------------

  matrix2(size_t m, size_t n)
  { resize(m,n); };

  matrix2(size_t m, size_t n, T D)
  { resize(m,n); for ( auto &i: m_data ) i = D; };

  matrix2(size_t m, size_t n, const T *D)
  { resize(m,n); for ( size_t i=0; i<size(); ++i ) m_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator matrix2<U> ()
  {
    matrix2<U> out(shape(0),shape(1));

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // resize matrix2
  // -------------

  void resize(size_t m, size_t n)
  {
    m_m = m;
    m_n = n;

    m_data.resize(m_m*m_n);
  }

  // reshape
  // -------

  void reshape(size_t m, size_t n)
  {
    assert( m_m*m_n == m*n );

    resize(m,n);
  }

  // operator[] : direct storage access
  // ----------------------------------

  T& operator[](size_t i)
  { return m_data[i]; };

  const T& operator[](size_t i) const
  { return m_data[i]; };

  // operator() : indices along each dimension
  // -----------------------------------------

  T& operator()(size_t a, size_t b)
  { return m_data[a*m_n+b]; };

  const T& operator()(size_t a, size_t b) const
  { return m_data[a*m_n+b]; };

  // arithmetic operators: matrix2 ?= matrix2
  // --------------------------------------

  matrix2<T>& operator*= (const matrix2<T> &B)
  {
    assert( shape(0) == B.shape(0) );
    assert( shape(1) == B.shape(1) );

    for ( size_t i = 0 ; i < m_m*m_n ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  matrix2<T>& operator/= (const matrix2<T> &B)
  {
    assert( shape(0) == B.shape(0) );
    assert( shape(1) == B.shape(1) );

    for ( size_t i = 0 ; i < m_m*m_n ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  matrix2<T>& operator+= (const matrix2<T> &B)
  {
    assert( shape(0) == B.shape(0) );
    assert( shape(1) == B.shape(1) );

    for ( size_t i = 0 ; i < m_m*m_n ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  matrix2<T>& operator-= (const matrix2<T> &B)
  {
    assert( shape(0) == B.shape(0) );
    assert( shape(1) == B.shape(1) );

    for ( size_t i = 0 ; i < m_m*m_n ; ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: matrix2 ?= scalar
  // --------------------------------------

  matrix2<T>& operator*= (T B)
  {
    for ( auto &i : m_data )
      i *= B;

    return *this;
  };

  matrix2<T>& operator/= (T B)
  {
    for ( auto &i : m_data )
      i /= B;

    return *this;
  };

  matrix2<T>& operator+= (T B)
  {
    for ( auto &i : m_data )
      i += B;

    return *this;
  };

  matrix2<T>& operator-= (T B)
  {
    for ( auto &i : m_data )
      i -= B;

    return *this;
  };

  // iterators / pointer
  // -------------------

  const T* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // return shape array [ndim]
  // -------------------------

  std::vector<size_t> shape() const
  {
    std::vector<size_t> ret(2);

    ret[0] = m_m;
    ret[1] = m_n;

    return ret;
  };

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return m_m;
    if ( i == 1 ) return m_n;

    assert( false );
  }

  // return strides array [ndim]
  // ---------------------------

  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> ret(2);

    ret[0] = m_n;
    ret[1] = 1;

    if ( bytes )
      for ( size_t i = 0 ; i < 2 ; ++i )
        ret[i] *= sizeof(T);

    return ret;
  };

  // return size
  // -----------

  size_t size() const { return m_m*m_n; };
  size_t ndim() const { return 2;       };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( auto &i : m_data )
      out += i;

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(m_m*m_n); };
  T      min () const { return *std::min_element(m_data.begin(),m_data.end()); };
  T      max () const { return *std::max_element(m_data.begin(),m_data.end()); };

  // initialize to zero/one/constant
  // -------------------------------

  void setConstant(T D) { for ( auto &i : m_data ) i = D;                 };
  void setZero    (   ) { for ( auto &i : m_data ) i = static_cast<T>(0); };
  void setOnes    (   ) { for ( auto &i : m_data ) i = static_cast<T>(0); };
  void zeros      (   ) { for ( auto &i : m_data ) i = static_cast<T>(0); };
  void ones       (   ) { for ( auto &i : m_data ) i = static_cast<T>(1); };

  // print to screen
  // ---------------

  void printf(std::string fmt) const
  {
    std::vector<size_t> s = strides();

    for ( size_t h = 0 ; h < shape(0) ; ++h )
    {
      for ( size_t i = 0 ; i < shape(1)-1 ; ++i )
        std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]]);

      std::printf((fmt+";\n").c_str(),m_data[h*s[0]+(shape()[1]-1)*s[1]]);
    }
  }

}; // class matrix2

// arithmetic operators: matrix2 = matrix2 ? matrix2
// ----------------------------------------------

template<class T> matrix2<T> operator* (const matrix2<T> &A, const matrix2<T> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template<class T> matrix2<T> operator/ (const matrix2<T> &A, const matrix2<T> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template<class T> matrix2<T> operator+ (const matrix2<T> &A, const matrix2<T> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template<class T> matrix2<T> operator- (const matrix2<T> &A, const matrix2<T> &B)
{
  assert( A.shape(0) == B.shape(0) );
  assert( A.shape(1) == B.shape(1) );

  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: matrix2 = matrix2 ? scalar
// ----------------------------------------------

template<class T> matrix2<T> operator* (const matrix2<T> &A, const T &B)
{
  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template<class T> matrix2<T> operator/ (const matrix2<T> &A, const T &B)
{
  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template<class T> matrix2<T> operator+ (const matrix2<T> &A, const T &B)
{
  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template<class T> matrix2<T> operator- (const matrix2<T> &A, const T &B)
{
  matrix2<T> C(A.shape(0),A.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: matrix2 = scalar ? matrix2
// ----------------------------------------------

template<class T> matrix2<T> operator* (const T &A, const matrix2<T> &B)
{
  matrix2<T> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template<class T> matrix2<T> operator/ (const T &A, const matrix2<T> &B)
{
  matrix2<T> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template<class T> matrix2<T> operator+ (const T &A, const matrix2<T> &B)
{
  matrix2<T> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template<class T> matrix2<T> operator- (const T &A, const matrix2<T> &B)
{
  matrix2<T> C(B.shape(0),B.shape(1));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class T>
std::ostream& operator<<(std::ostream& out, matrix2<T>& src)
{
  for ( size_t i = 0 ; i < src.shape(0) ; ++i )
  {
    for ( size_t j = 0 ; j < src.shape(1)-1 ; ++j )
      out << src(i,j) << ", ";

    out << src(i,src.shape()[1]-1) << "; " << std::endl;
  }

  return out;
}

}} // namespace cppmat::periodic

#endif

