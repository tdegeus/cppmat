/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX2_H
#define CPPMAT_TINY_MATRIX2_H

#include "macros.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::matrix2
// =================================================================================================

template <class X, size_t m, size_t n> class matrix2
{
private:

  X      m_container[m*n];  // data container
  X     *m_data;            // pointer to container (may point outside)
  size_t m_size=m*n;        // total number of entries

public:

  // constructors
  matrix2();
  matrix2(X D);
  matrix2(const X *D);
  matrix2(const matrix2<X,m,n> &D);

  // copy constructor
  matrix2<X,m,n>& operator= (const matrix2<X,m,n> &D);

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
    typename U,typename V=X, size_t M, size_t N,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator matrix2<U,M,N> ()
  {
    matrix2<U,M,N> out;

    for ( size_t i = 0 ; i < M*N ; ++i )
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

  X& operator()(size_t a, size_t b)
  { return m_data[a*n+b]; }

  const X& operator()(size_t a, size_t b) const
  { return m_data[a*n+b]; }

  // arithmetic operators: matrix2 ?= matrix2
  // --------------------------------------

  matrix2<X,m,n>& operator*= (const matrix2<X,m,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  matrix2<X,m,n>& operator/= (const matrix2<X,m,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  matrix2<X,m,n>& operator+= (const matrix2<X,m,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  matrix2<X,m,n>& operator-= (const matrix2<X,m,n> &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: matrix2 ?= scalar
  // --------------------------------------

  matrix2<X,m,n>& operator*= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  matrix2<X,m,n>& operator/= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  matrix2<X,m,n>& operator+= (X B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  matrix2<X,m,n>& operator-= (X B)
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
    std::vector<size_t> ret(2);

    ret[0] = m;
    ret[1] = n;

    return ret;
  }

  // return shape in one direction
  // -----------------------------

  size_t shape(size_t i) const
  {
    if ( i == 0 ) return m;
    if ( i == 1 ) return n;

    assert( false );
    return 0;
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
        ret[i] *= sizeof(X);

    return ret;
  }

  // return size
  // -----------

  size_t size() const { return m_size; }
  size_t ndim() const { return 2;      }

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

    for ( size_t h = 0 ; h < shape(0) ; ++h )
    {
      for ( size_t i = 0 ; i < shape(1)-1 ; ++i )
        std::printf((fmt+",").c_str(),m_data[h*s[0]+i*s[1]]);

      std::printf((fmt+";\n").c_str(),m_data[h*s[0]+(shape(1)-1)*s[1]]);
    }
  }

}; // class matrix2

// arithmetic operators: matrix2 = matrix2 ? matrix2
// -------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: matrix2 = matrix2 ? scalar
// ------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const X &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: matrix2 = scalar ? matrix2
// ------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator* (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator/ (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator+ (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template <class X, size_t m, size_t n>
matrix2<X,m,n> operator- (const X &A, const matrix2<X,m,n> &B)
{
  matrix2<X,m,n> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// print to "std::cout"
// --------------------

template <class X, size_t m, size_t n>
std::ostream& operator<<(std::ostream& out, matrix2<X,m,n>& src)
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
// constructors
// =================================================================================================

template <class X, size_t m, size_t n>
matrix2<X,m,n>::matrix2()
{
  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n>::matrix2(X D)
{
  // copy input
  for ( size_t i = 0; i < m_size ; ++i ) m_container[i] = D;

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n>::matrix2(const X *D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// -------------------------------------------------------------------------------------------------

template <class X, size_t m, size_t n>
matrix2<X,m,n>::matrix2(const matrix2<X,m,n> &D)
{
  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i ) m_container[i] = D[i];

  // point to local data container
  m_data = &m_container[0];
}

// =================================================================================================
// copy constructor
// =================================================================================================

template <class X, size_t m, size_t n>
matrix2<X,m,n>& matrix2<X,m,n>::operator= (const matrix2<X,m,n> &D)
{
  // - copy input
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_container[i] = D[i];
  // - point to local data container
  m_data = &m_container[0];
  // - return pointer to current instance
  return *this;
}

// =================================================================================================

} // namespace tiny
} // namespace cppmat

#endif

