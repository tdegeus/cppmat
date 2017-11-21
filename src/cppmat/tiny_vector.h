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

template <class X, size_t n> class vector
{
private:

  X      m_container[n];    // data container
  X     *m_data;            // pointer to container (may point outside)
  size_t m_size=n;          // total number of entries

public:

  // constructors
  // ------------

  vector()
  {
    // - point to local data container
    m_data = &m_container[0];
  }

  vector(X D)
  {
    // - copy input
    for ( size_t i = 0; i < m_size ; ++i )
      m_container[i] = D;
    // - point to local data container
    m_data = &m_container[0];
  }

  vector(const X *D)
  {
    // - copy input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_container[i] = D[i];
    // - point to local data container
    m_data = &m_container[0];
  }

  vector(const vector<X,n> &D)
  {
    // - copy input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_container[i] = D[i];
    // - point to local data container
    m_data = &m_container[0];
  }

  vector<X,n>& operator= (const vector<X,n> &D)
  {
    // - copy input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_container[i] = D[i];
    // - point to local data container
    m_data = &m_container[0];
    // - return pointer to current instance
    return *this;
  }

  // map external pointer
  // --------------------

  // raw pointer
  void map(X *D)
  {
    m_data = D;
  }

  // constructor to copy + change data type
  // --------------------------------------

  template<\
    typename U,typename V=X, size_t N,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator vector<U,N> ()
  {
    vector<U,N> out;

    for ( size_t i = 0 ; i < N ; ++i )
      out[i] = static_cast<U>(m_data[i]);

    return out;
  }

  // operator[] : direct storage access
  // ----------------------------------

  X& operator[](size_t i)
  { return m_data[i]; }

  const X& operator[](size_t i) const
  { return m_data[i]; }

  // operator() : indices along each dimension
  // -----------------------------------------

  X& operator()(size_t a)
  { return m_data[a]; }

  const X& operator()(size_t a) const
  { return m_data[a]; }

  // arithmetic operators: vector ?= vector
  // --------------------------------------

  vector<X,n>& operator*= (const vector<X,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  vector<X,n>& operator/= (const vector<X,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  vector<X,n>& operator+= (const vector<X,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  vector<X,n>& operator-= (const vector<X,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<X,n>& operator*= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  vector<X,n>& operator/= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  vector<X,n>& operator+= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  vector<X,n>& operator-= (X B)
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

    ret[0] = n;

    return ret;
  }

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return n;

    assert( false );
    return 0;
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

template <class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: vector = vector ? scalar
// ----------------------------------------------

template <class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const X &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
// ----------------------------------------------

template <class X, size_t n>
inline vector<X,n> operator* (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator/ (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator+ (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template <class X, size_t n>
inline vector<X,n> operator- (const X &A, const vector<X,n> &B)
{
  vector<X,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class X, size_t n>
inline std::ostream& operator<<(std::ostream& out, vector<X,n>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

} // namespace tiny
} // namespace cppmat

#endif

