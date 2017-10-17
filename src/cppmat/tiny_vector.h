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

template <class T, size_t n> class vector
{
private:

  T m_data[n]; // data array

public:

  // (copy) constructor
  // ------------------

  vector                 (const vector<T,n> &) = default;
  vector<T,n>& operator= (const vector<T,n> &) = default;
  vector<T,n>(){};

  // explicit constructors
  // ---------------------

  vector(T D)
  { for ( size_t i = 0; i < n ; ++i ) m_data[i] = D; };

  vector(const T *D)
  { for ( size_t i = 0; i < n ; ++i ) m_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T, size_t N,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator vector<U,N> ()
  {
    vector<U,N> out;

    for ( size_t i = 0 ; i < N ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
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

  vector<T,n>& operator*= (const vector<T,n> &B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  vector<T,n>& operator/= (const vector<T,n> &B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  vector<T,n>& operator+= (const vector<T,n> &B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  vector<T,n>& operator-= (const vector<T,n> &B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<T,n>& operator*= (T B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] *= B;

    return *this;
  };

  vector<T,n>& operator/= (T B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] /= B;

    return *this;
  };

  vector<T,n>& operator+= (T B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] += B;

    return *this;
  };

  vector<T,n>& operator-= (T B)
  {
    for ( size_t i = 0 ; i < n ; ++i )
      m_data[i] -= B;

    return *this;
  };

  // iterators / pointer
  // -------------------

  const T* data () const { return &m_data[0]; };
  auto     begin()       { return &m_data[0]; };
  auto     end  ()       { return &m_data[n-1]; };

  // return shape array [ndim]
  // -------------------------

  std::vector<size_t> shape() const
  {
    std::vector<size_t> ret(1);

    ret[0] = n;

    return ret;
  };

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return n;

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

  size_t size() const { return n; };
  size_t ndim() const { return 1; };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( size_t i = 0 ; i < n ; ++i )
      out += m_data[i];

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(n); };
  T      min () const { return *std::min_element(begin(),end()); };
  T      max () const { return *std::max_element(begin(),end()); };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = static_cast<T>(0); };
  void ones () { for ( size_t i = 0 ; i < n ; ++i ) m_data[i] = static_cast<T>(1); };

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

template <class T, size_t n>
vector<T,n> operator* (const vector<T,n> &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator/ (const vector<T,n> &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator+ (const vector<T,n> &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator- (const vector<T,n> &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: vector = vector ? scalar
// ----------------------------------------------

template <class T, size_t n>
vector<T,n> operator* (const vector<T,n> &A, const T &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template <class T, size_t n>
vector<T,n> operator/ (const vector<T,n> &A, const T &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template <class T, size_t n>
vector<T,n> operator+ (const vector<T,n> &A, const T &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template <class T, size_t n>
vector<T,n> operator- (const vector<T,n> &A, const T &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
// ----------------------------------------------

template <class T, size_t n>
vector<T,n> operator* (const T &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator/ (const T &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator+ (const T &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template <class T, size_t n>
vector<T,n> operator- (const T &A, const vector<T,n> &B)
{
  vector<T,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class T, size_t n>
std::ostream& operator<<(std::ostream& out, vector<T,n>& src)
{
  for ( size_t i = 0 ; i < src.shape(0)-1 ; ++i )
    out << src(i) << " , ";

  out << src(src.shape(0)-1) << std::endl;

  return out;
}

} // namespace tiny
} // namespace cppmat

#endif

