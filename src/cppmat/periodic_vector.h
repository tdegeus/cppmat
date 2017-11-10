/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_VECTOR_H
#define CPPMAT_PERIODIC_VECTOR_H

#include "macros.h"

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template <class X> class vector
{
private:

  std::vector<X> m_container;   // data container
  X             *m_data;        // pointer to data container (may point outside)
  size_t         m_n=0;         // number of columns
  size_t         m_size=0;      // total size
  bool           m_owner=true;  // signal if m_data pointer to "m_container"

public:

  // constructors
  // ------------

  vector(){}

  vector(size_t n){ resize(n); }

  vector(size_t n, X D)
  {
    resize(n);

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D;
  }

  vector(size_t n, const X *D)
  {
    resize(n);

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t n, X *D)
  {
    // - release 'ownership' of data
    m_owner = false;

    // - change size settings without expanding "m_container"
    resize(n);

    // - point to input pointer
    m_data = D;
  }

  // copy from external data array
  // -----------------------------

  void copy(size_t n, const X *D)
  {
    assert( m_owner );

    resize(n);

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // constructor to copy + change data type
  // --------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator vector<U> ()
  {
    vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> ()
  {
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // resize vector
  // -------------

  void resize(size_t n)
  {
    m_n    = n;
    m_size = n;

    if ( m_owner )
    {
      m_container.resize(m_size);
      m_data = &m_container[0];
    }
  }

  // operator[] : direct storage access
  // ----------------------------------

  X& operator[](size_t i)
  { return m_data[i]; }

  const X& operator[](size_t i) const
  { return m_data[i]; };

  // operator() : indices along each dimension
  // -----------------------------------------

  X& operator()(size_t a)
  { return m_data[a]; };

  const X& operator()(size_t a) const
  { return m_data[a]; };

  // operator() : indices along each dimension, allowing for "i < 0" or "i >= N"
  // ---------------------------------------------------------------------------

  X& operator()(int a)
  {
    a = ( a < 0 ) ? a + static_cast<int>(m_size) : ( a >= static_cast<int>(m_size) ) ? a - static_cast<int>(m_size) : a ;

    return m_data[a];
  }

  const X& operator()(int a) const
  {
    a = ( a < 0 ) ? a + static_cast<int>(m_size) : ( a >= static_cast<int>(m_size) ) ? a - static_cast<int>(m_size) : a ;

    return m_data[a];
  }

  // arithmetic operators: vector ?= vector
  // --------------------------------------

  vector<X>& operator*= (const vector<X> &B)
  {
    assert( m_n == B.shape(0) );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  vector<X>& operator/= (const vector<X> &B)
  {
    assert( m_n == B.shape(0) );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  vector<X>& operator+= (const vector<X> &B)
  {
    assert( m_n == B.shape(0) );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  vector<X>& operator-= (const vector<X> &B)
  {
    assert( m_n == B.shape(0) );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<X>& operator*= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  vector<X>& operator/= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  vector<X>& operator+= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  vector<X>& operator-= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B;

    return *this;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // return shape array [ndim]
  // -------------------------

  std::vector<size_t> shape() const
  {
    std::vector<size_t> ret(1);

    ret[0] = m_n;

    return ret;
  }

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
      ret[0] *= sizeof(X);

    return ret;
  }

  // return size
  // -----------

  size_t size() const { return m_size; }
  size_t ndim() const { return 1;      }

  // minimum / maximum / mean / sum
  // ------------------------------

  X sum() const
  {
    X out = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out += m_data[i];

    return out;
  }

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(m_size); }
  X      min () const { return *std::min_element(begin(),end()); }
  X      max () const { return *std::max_element(begin(),end()); }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

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

template<class X> vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template<class X> vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template<class X> vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template<class X> vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: vector = vector ? scalar
// ----------------------------------------------

template<class X> vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template<class X> vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template<class X> vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template<class X> vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
// ----------------------------------------------

template<class X> vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template<class X> vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template<class X> vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template<class X> vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

}} // namespace cppmat::periodic

#endif

