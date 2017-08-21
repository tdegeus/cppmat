/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

namespace cppmat {

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template <class T> class matrix
{
private:

  std::vector<T>      m_data;    // data array
  std::vector<size_t> m_shape;   // number of entries in each dimensions
  std::vector<size_t> m_strides; // stride length for each index
  // Note:
  // - Number of dimensions determined from "m_shape.size()".
  // - "m_strides.size()" can be bigger than "m_shape.size()" to allow reading with more indices.
  //   E.g. after using "atleast_3d" also a 2D matrix can be read using 3 indices.

public:

  // (copy) constructor
  // ------------------

  matrix               (const matrix<T> &) = default;
  matrix<T>& operator= (const matrix<T> &) = default;
  matrix<T>(){};

  // explicit constructors
  // ---------------------

  matrix(const std::vector<size_t> &shape)
  { resize(shape); };

  matrix(const std::vector<size_t> &shape, T D)
  { resize(shape); for ( auto &i: m_data ) i = D; };

  matrix(const std::vector<size_t> &shape, const T *D)
  { resize(shape); for ( size_t i=0; i<size(); ++i ) m_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator matrix<U> ()
  {
    matrix<U> out(shape());

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // resize matrix
  // -------------

  void resize(const std::vector<size_t> &shape)
  {
    if ( shape.size()<1 )
      throw std::runtime_error("Input should be >= 1-D");

    m_shape  .resize(shape.size());
    m_strides.resize(shape.size());
    size_t n = shape[0];

    for ( size_t i=1; i<shape.size(); ++i ) { n           *= shape[i]; }
    for ( size_t i=0; i<shape.size(); ++i ) { m_shape  [i] = shape[i]; }
    for ( size_t i=0; i<shape.size(); ++i ) { m_strides[i] = 1;        }

    for ( size_t i=0; i<shape.size(); ++i )
      for ( size_t j=i+1; j<shape.size(); ++j )
        m_strides[i] *= m_shape[j];

    m_data.resize(n);
  }

  // convert matrix to view it in d >= nd
  // ------------------------------------

  void atleast_1d()         { while ( m_strides.size()<1 ) { m_strides.push_back(1); } }
  void atleast_2d()         { while ( m_strides.size()<2 ) { m_strides.push_back(1); } }
  void atleast_3d()         { while ( m_strides.size()<3 ) { m_strides.push_back(1); } }
  void atleast_nd(size_t n) { while ( m_strides.size()<n ) { m_strides.push_back(1); } }

  // index operators (now up to 6-d, extend if needed)
  // -------------------------------------------------

  T& operator[](size_t i)
  { return m_data[i]; };

  T& operator()(size_t a)
  { return m_data[a*m_strides[0]]; };

  T& operator()(size_t a, size_t b)
  { return m_data[a*m_strides[0]+b*m_strides[1]]; };

  T& operator()(size_t a, size_t b, size_t c)
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d)
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]]; };

  T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]]; };

  const T& operator[](size_t i) const
  { return m_data[i]; };

  const T& operator()(size_t a) const
  { return m_data[a*m_strides[0]]; };

  const T& operator()(size_t a, size_t b) const
  { return m_data[a*m_strides[0]+b*m_strides[1]]; };

  const T& operator()(size_t a, size_t b, size_t c) const
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]]; };

  const T& operator()(size_t a, size_t b, size_t c, size_t d) const
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]]; };

  const T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]]; };

  const T& operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
  { return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]]; };

  // arithmetic operators: matrix ?= matrix
  // --------------------------------------

  matrix<T>& operator*= (const matrix<T> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < size() ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  matrix<T>& operator/= (const matrix<T> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < size() ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  matrix<T>& operator+= (const matrix<T> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < size() ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  matrix<T>& operator-= (const matrix<T> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < size() ; ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: matrix ?= scalar
  // --------------------------------------

  matrix<T>& operator*= (T B)
  {
    for ( auto &i : m_data )
      i *= B;

    return *this;
  };

  matrix<T>& operator/= (T B)
  {
    for ( auto &i : m_data )
      i /= B;

    return *this;
  };

  matrix<T>& operator+= (T B)
  {
    for ( auto &i : m_data )
      i += B;

    return *this;
  };

  matrix<T>& operator-= (T B)
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

  std::vector<size_t> shape(size_t nd=0) const
  {
    if ( nd==0 )
      nd = ndim();

    std::vector<size_t> ret(nd,1);

    for ( size_t i=0; i<nd; ++i )
      ret[i] = m_shape[i];

    return ret;
  };

  // return strides array [ndim]
  // ---------------------------

  std::vector<size_t> strides(bool bytes=false) const
  {
    size_t nd = ndim();
    std::vector<size_t> ret(nd);

    for ( size_t i=0; i<nd; ++i )
      ret[i] = m_strides[i];

    if ( bytes )
      for ( size_t i=0; i<nd; ++i )
        ret[i] *= sizeof(T);

    return ret;
  };

  // return size
  // -----------

  size_t size() const { return m_data .size(); };
  size_t ndim() const { return m_shape.size(); };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( auto i: m_data )
      out += i;

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(size()); };
  T      min () const { return *std::min_element(m_data.begin(),m_data.end()); };
  T      max () const { return *std::max_element(m_data.begin(),m_data.end()); };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( auto &i: m_data ) i = static_cast<T>(0); };
  void ones () { for ( auto &i: m_data ) i = static_cast<T>(1); };

  // print to screen
  // ---------------

  void printf(std::string fmt) const
  {
    std::vector<size_t> s = m_strides;
    if ( ndim()==1 ) {
      for ( size_t h=0; h<shape()[0]-1; ++h )
        std::printf((fmt+",").c_str(),m_data[h]);
      std::printf((fmt+"\n").c_str(),m_data[shape()[0]-1]);
    }
    else if ( ndim()==2 ) {
      for ( size_t h=0; h<shape()[0]; ++h ) {
        for ( size_t i=0; i<shape()[1]-1; ++i )
          std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]]);
        std::printf((fmt+";\n").c_str(),m_data[h*s[0]+(shape()[1]-1)*s[1]]);
      }
    }
    else if ( ndim()==3 ) {
      for ( size_t h=0; h<shape()[0]; ++h ) {
        for ( size_t i=0; i<shape()[1]; ++i ) {
          for ( size_t j=0; j<shape()[2]-1; ++j )
            std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]+j*s[2]]);
          std::printf((fmt+";\n").c_str(),m_data[h*s[0]+i*s[1]+(shape()[2]-1)*s[2]]);
        }
        if ( h<shape()[0]-1 )
          std::printf("\n");
      }
    }
  }

}; // class matrix

// arithmetic operators: matrix = matrix ? matrix
// ----------------------------------------------

template<class T> matrix<T> operator* (const matrix<T> &A, const matrix<T> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template<class T> matrix<T> operator/ (const matrix<T> &A, const matrix<T> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B[i];

  return C;
}

template<class T> matrix<T> operator+ (const matrix<T> &A, const matrix<T> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template<class T> matrix<T> operator- (const matrix<T> &A, const matrix<T> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: matrix = matrix ? scalar
// ----------------------------------------------

template<class T> matrix<T> operator* (const matrix<T> &A, const T &B)
{
  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B;

  return C;
}

template<class T> matrix<T> operator/ (const matrix<T> &A, const T &B)
{
  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B;

  return C;
}

template<class T> matrix<T> operator+ (const matrix<T> &A, const T &B)
{
  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B;

  return C;
}

template<class T> matrix<T> operator- (const matrix<T> &A, const T &B)
{
  matrix<T> C(A.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: matrix = scalar ? matrix
// ----------------------------------------------

template<class T> matrix<T> operator* (const T &A, const matrix<T> &B)
{
  matrix<T> C(B.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A * B[i];

  return C;
}

template<class T> matrix<T> operator/ (const T &A, const matrix<T> &B)
{
  matrix<T> C(B.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A / B[i];

  return C;
}

template<class T> matrix<T> operator+ (const T &A, const matrix<T> &B)
{
  matrix<T> C(B.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A + B[i];

  return C;
}

template<class T> matrix<T> operator- (const T &A, const matrix<T> &B)
{
  matrix<T> C(B.shape());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class T>
std::ostream& operator<<(std::ostream& out, matrix<T>& src)
{
  if ( src.ndim()==1 ) {
    for ( size_t i=0; i<src.shape()[0]-1; ++i )
      out << src(i) << " , ";
    out << src(src.shape()[0]-1) << std::endl;
  }
  else if ( src.ndim()==2 ) {
    for ( size_t i=0; i<src.shape()[0]; ++i ) {
      for ( size_t j=0; j<src.shape()[1]-1; ++j ) {
        out << src(i,j) << ", ";
      }
      out << src(i,src.shape()[1]-1) << "; " << std::endl;
    }
  }
  else if ( src.ndim()==3 ) {
    for ( size_t h=0; h<src.shape()[0]; ++h ) {
      for ( size_t i=0; i<src.shape()[1]; ++i ) {
        for ( size_t j=0; j<src.shape()[2]-1; ++j ) {
          out << src(h,i,j) << ", ";
        }
        out << src(i,src.shape()[1]-1) << "; " << std::endl;
      }
      if ( h<src.shape()[0]-1 )
        out << std::endl;
    }
  }

  return out;
}

} // namespace cppmat

#endif

