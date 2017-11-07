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

template <class T> class vector
{
private:

  std::vector<T> m_data;   // data array
  size_t         m_n=0;    // number of columns

public:

  // (copy) constructor
  // ------------------

  vector               (const vector<T> &) = default;
  vector<T>& operator= (const vector<T> &) = default;
  vector<T>(){};

  // explicit constructors
  // ---------------------

  vector(size_t n)
  { resize(n); };

  vector(size_t n, T D)
  { resize(n); for ( auto &i: m_data ) i = D; };

  vector(size_t n, const T *D)
  { resize(n); for ( size_t i=0; i<size(); ++i ) m_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator vector<U> ()
  {
    vector<U> out(shape(0));

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // resize vector
  // -------------

  void resize(size_t n)
  {
    m_n = n;

    m_data.resize(m_n);
  }

  // operator[] : direct storage access
  // ----------------------------------

  T& operator[](size_t i)
  { return m_data[i]; };

  const T& operator[](size_t i) const
  { return m_data[i]; };

  // operator() : indices along each dimension
  // -----------------------------------------

  T& operator()(size_t a)
  { return m_data[a]; };

  const T& operator()(size_t a) const
  { return m_data[a]; };

  // arithmetic operators: vector ?= vector
  // --------------------------------------

  vector<T>& operator*= (const vector<T> &B)
  {
    assert( shape(0) == B.shape(0) );

    for ( size_t i = 0 ; i < m_n ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  vector<T>& operator/= (const vector<T> &B)
  {
    assert( shape(0) == B.shape(0) );

    for ( size_t i = 0 ; i < m_n ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  vector<T>& operator+= (const vector<T> &B)
  {
    assert( shape(0) == B.shape(0) );

    for ( size_t i = 0 ; i < m_n ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  vector<T>& operator-= (const vector<T> &B)
  {
    assert( shape(0) == B.shape(0) );

    for ( size_t i = 0 ; i < m_n ; ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<T>& operator*= (T B)
  {
    for ( auto &i : m_data )
      i *= B;

    return *this;
  };

  vector<T>& operator/= (T B)
  {
    for ( auto &i : m_data )
      i /= B;

    return *this;
  };

  vector<T>& operator+= (T B)
  {
    for ( auto &i : m_data )
      i += B;

    return *this;
  };

  vector<T>& operator-= (T B)
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
    std::vector<size_t> ret(1);

    ret[0] = m_n;

    return ret;
  };

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return m_n;

    assert( false );
  }

  // return strides array [ndim]
  // ---------------------------

  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> ret(1);

    ret[0] = 1;

    if ( bytes )
      ret[0] *= sizeof(T);

    return ret;
  };

  // return size
  // -----------

  size_t size() const { return m_n; };
  size_t ndim() const { return 2;   };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( auto &i : m_data )
      out += i;

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(m_n); };
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

    for ( size_t h = 0 ; h < shape(0)-1 ; ++h )
      std::printf((fmt+",").c_str(),m_data[h]);

    std::printf((fmt+"\n").c_str(),m_data[shape(0)-1]);
  }

}; // class vector

// arithmetic operators: vector = vector ? vector
// ----------------------------------------------

template<class T> vector<T> operator* (const vector<T> &A, const vector<T> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template<class T> vector<T> operator/ (const vector<T> &A, const vector<T> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template<class T> vector<T> operator+ (const vector<T> &A, const vector<T> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template<class T> vector<T> operator- (const vector<T> &A, const vector<T> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: vector = vector ? scalar
// ----------------------------------------------

template<class T> vector<T> operator* (const vector<T> &A, const T &B)
{
  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template<class T> vector<T> operator/ (const vector<T> &A, const T &B)
{
  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template<class T> vector<T> operator+ (const vector<T> &A, const T &B)
{
  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template<class T> vector<T> operator- (const vector<T> &A, const T &B)
{
  vector<T> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
// ----------------------------------------------

template<class T> vector<T> operator* (const T &A, const vector<T> &B)
{
  vector<T> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template<class T> vector<T> operator/ (const T &A, const vector<T> &B)
{
  vector<T> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template<class T> vector<T> operator+ (const T &A, const vector<T> &B)
{
  vector<T> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template<class T> vector<T> operator- (const T &A, const vector<T> &B)
{
  vector<T> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class T>
std::ostream& operator<<(std::ostream& out, vector<T>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

} // namespace cppmat

#endif

