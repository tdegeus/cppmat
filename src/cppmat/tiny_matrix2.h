/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX2_H
#define CPPMAT_TINY_MATRIX2_H

#include "macros.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::matrix2
// =================================================================================================

template <class T, size_t m, size_t n> class matrix2
{
private:

  T m_data[m*n]; // data array

public:

  // (copy) constructor
  // ------------------

  matrix2                   (const matrix2<T,m,n> &) = default;
  matrix2<T,m,n>& operator= (const matrix2<T,m,n> &) = default;
  matrix2<T,m,n>(){};

  // explicit constructors
  // ---------------------

  matrix2(T D)
  { for ( size_t i = 0; i < m*n ; ++i ) m_data[i] = D; };

  matrix2(const T *D)
  { for ( size_t i = 0; i < m*n ; ++i ) m_data[i] = D[i]; };

  // constructor to copy + change data type
  // --------------------------------------

  template<typename U,typename V=T, size_t M, size_t N,\
    typename=typename std::enable_if<std::is_convertible<T,U>::value>::type>
  operator matrix2<U,M,N> ()
  {
    matrix2<U,M,N> out;

    for ( size_t i = 0 ; i < M*N ; ++i )
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

  T& operator()(size_t a, size_t b)
  { return m_data[a*n+b]; };

  const T& operator()(size_t a, size_t b) const
  { return m_data[a*n+b]; };

  // arithmetic operators: matrix2 ?= matrix2
  // --------------------------------------

  matrix2<T,m,n>& operator*= (const matrix2<T,m,n> &B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] *= B[i];

    return *this;
  };

  matrix2<T,m,n>& operator/= (const matrix2<T,m,n> &B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] /= B[i];

    return *this;
  };

  matrix2<T,m,n>& operator+= (const matrix2<T,m,n> &B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] += B[i];

    return *this;
  };

  matrix2<T,m,n>& operator-= (const matrix2<T,m,n> &B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: matrix2 ?= scalar
  // --------------------------------------

  matrix2<T,m,n>& operator*= (T B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] *= B;

    return *this;
  };

  matrix2<T,m,n>& operator/= (T B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] /= B;

    return *this;
  };

  matrix2<T,m,n>& operator+= (T B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] += B;

    return *this;
  };

  matrix2<T,m,n>& operator-= (T B)
  {
    for ( size_t i = 0 ; i < m*n ; ++i )
      m_data[i] -= B;

    return *this;
  };

  // iterators / pointer
  // -------------------

  const T* data () const { return &m_data[0]; };
  auto     begin()       { return &m_data[0]; };
  auto     end  ()       { return &m_data[m*n-1]; };

  // return shape array [ndim]
  // -------------------------

  std::vector<size_t> shape() const
  {
    std::vector<size_t> ret(2);

    ret[0] = m;
    ret[1] = n;

    return ret;
  };

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return m;
    if ( i == 1 ) return n;

    assert( false );
  }

  // return strides array [ndim]
  // ---------------------------

  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> ret(2);

    ret[0] = n;
    ret[1] = 1;

    if ( bytes )
      for ( size_t i = 0 ; i < 2 ; ++i )
        ret[i] *= sizeof(T);

    return ret;
  };

  // return size
  // -----------

  size_t size() const { return m*n; };
  size_t ndim() const { return 2;   };

  // minimum / maximum / mean / sum
  // ------------------------------

  T sum() const
  {
    T out = static_cast<T>(0);

    for ( size_t i = 0 ; i < m*n ; ++i )
      out += m_data[i];

    return out;
  };

  double mean() const { return static_cast<double>(this->sum())/static_cast<double>(m*n); };
  T      min () const { return *std::min_element(begin(),end()); };
  T      max () const { return *std::max_element(begin(),end()); };

  // initialize to zero/one
  // ----------------------

  void zeros() { for ( size_t i = 0 ; i < m*n ; ++i ) m_data[i] = static_cast<T>(0); };
  void ones () { for ( size_t i = 0 ; i < m*n ; ++i ) m_data[i] = static_cast<T>(1); };

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

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator* (const matrix2<T,m,n> &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator/ (const matrix2<T,m,n> &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator+ (const matrix2<T,m,n> &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator- (const matrix2<T,m,n> &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: matrix2 = matrix2 ? scalar
// ----------------------------------------------

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator* (const matrix2<T,m,n> &A, const T &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator/ (const matrix2<T,m,n> &A, const T &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator+ (const matrix2<T,m,n> &A, const T &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator- (const matrix2<T,m,n> &A, const T &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: matrix2 = scalar ? matrix2
// ----------------------------------------------

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator* (const T &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator/ (const T &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator+ (const T &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template<class T, size_t m, size_t n>
matrix2<T,m,n> operator- (const T &A, const matrix2<T,m,n> &B)
{
  matrix2<T,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class T, size_t m, size_t n>
std::ostream& operator<<(std::ostream& out, matrix2<T,m,n>& src)
{
  for ( size_t i = 0 ; i < src.shape(0) ; ++i )
  {
    for ( size_t j = 0 ; j < src.shape(1)-1 ; ++j )
      out << src(i,j) << ", ";

    out << src(i,src.shape()[1]-1) << "; " << std::endl;
  }

  return out;
}

} // namespace tiny
} // namespace cppmat

#endif

