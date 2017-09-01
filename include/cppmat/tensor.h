/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef TENSOR_H
#define TENSOR_H

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

namespace cppmat {

// =================================================================================================
// forward declaration
// =================================================================================================

template<class X> class tensor4 ;
template<class X> class tensor2 ;
template<class X> class tensor2s;
template<class X> class tensor2d;
template<class X> class vector  ;

// =================================================================================================
// return strides (generic routine used by all tensor-classes)
// - defines how much to skip per index, e.g. for a tensor rank 2 of dimension 3: [3,1]
// - if bytes == true the definition is is bytes
// =================================================================================================

template<class X> std::vector<size_t> inline _strides(size_t rank, size_t ndim, bool bytes=false)
{
  std::vector<size_t> out(rank,1);

  for ( size_t i=0; i<rank; ++i )
    for ( size_t j=i+1; j<rank; ++j )
      out[i] *= ndim;

  if ( bytes )
    for ( auto &i: out )
      i *= sizeof(X);

  return out;
}

// =================================================================================================
// cppmat::tensor4
// =================================================================================================

template<class X> class tensor4
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd;   // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor4(){};

  // explicit constructor: set correct size (WARNING: data not initialized)
  tensor4(size_t nd ) { resize(nd); };

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor4(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D;    };
  tensor4(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D[i]; };

  // change number of dimensions (WARNING: data not initialized)
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd*nd*nd*nd); };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(4,m_nd,bytes); };

  // copy constructor
  // ----------------

  // copy "tensor4" -> "tensor4" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor4<U> () const
  {
    tensor4<U> out(m_nd);

    for ( size_t i=0; i<size(); ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // norm
  // ----

  X norm() const { X C = static_cast<X>(0); for ( auto &i : m_data ) C += std::abs(i); return C; }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( auto &i : m_data ) i = static_cast<X>(0); };
  void ones        (     ) { for ( auto &i : m_data ) i = static_cast<X>(1); };
  void setConstant ( X D ) { for ( auto &i : m_data ) i = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor4 <X> inline ddot(const tensor4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2 <X> inline ddot(const tensor2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2 <X> inline ddot(const tensor2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2 <X> inline ddot(const tensor2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor4 <X> inline T   (                    ) const; // transposition   : B_lkji = A_ijkl
  tensor4 <X> inline RT  (                    ) const; // transposition   : B_ijlk = A_ijkl
  tensor4 <X> inline LT  (                    ) const; // transposition   : B_jikl = A_ijkl

  // index operators
  // ---------------

  X& operator[](size_t i)
  { return m_data[i]; };

  const X& operator[](size_t i) const
  { return m_data[i]; };

  X& operator()(size_t i, size_t j, size_t k, size_t l)
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; };

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; };

  // arithmetic operators: tensor4 ?= tensor4
  // ----------------------------------------

  tensor4<X>& operator*= (const tensor4<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] *= B[i];

    return *this;
  };

  tensor4<X>& operator/= (const tensor4<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] /= B[i];

    return *this;
  };

  tensor4<X>& operator+= (const tensor4<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] += B[i];

    return *this;
  };

  tensor4<X>& operator-= (const tensor4<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor4 ?= tensor4
  // ----------------------------------------

  tensor4<X>& operator*= (const X &B)
  {
    for ( auto &i: m_data )
      i *= B;

    return *this;
  };

  tensor4<X>& operator/= (const X &B)
  {
    for ( auto &i: m_data )
      i /= B;

    return *this;
  };

  tensor4<X>& operator+= (const X &B)
  {
    for ( auto &i: m_data )
      i += B;

    return *this;
  };

  tensor4<X>& operator-= (const X &B)
  {
    for ( auto &i: m_data )
      i -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor4<X> &B )
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i<size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class tensor4

// arithmetic operators: tensor4 = tensor4 ? tensor4
// -------------------------------------------------

template <class X> tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X> tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: tensor4 = tensor4 ? scalar
// ------------------------------------------------

template <class X> tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B;

  return C; }

template <class X> tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B;

  return C; }

template <class X> tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B;

  return C; }

template <class X> tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B;

  return C; }

// arithmetic operators: tensor4 = scalar ? tensor4
// ------------------------------------------------

template <class X> tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A * B[i];

  return C;
}

template <class X> tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A / B[i];

  return C;
}

template <class X> tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A + B[i];

  return C;
}

template <class X> tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// cppmat::tensor2
// =================================================================================================

template<class X> class tensor2
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd;   // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2(){};

  // explicit constructor: set correct size (WARNING: data not initialized)
  tensor2(size_t nd ) { resize(nd); };

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor2(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D;    };
  tensor2(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D[i]; };

  // change number of dimensions (WARNING: data not initialized)
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd*nd); };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,m_nd,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2" -> "tensor2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd);

    for ( size_t i=0; i<size(); ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // convert "tensor2 -> tensor2s"
  // WARNING: the output is symmetrized: "out(i,j) = ( this(i,j) + this(j,i) ) / 2."
  tensor2s<X> astensor2s()
  {
    tensor2s<X> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i )
      for ( size_t j=i; j<m_nd; ++j )
        out[ i*m_nd - (i-1)*i/2 + j - i ] =
        ( m_data[ i*m_nd + j ] + m_data[ j*m_nd + i ] ) / static_cast<X>(2);

    return out;
  }

  // convert "tensor2 -> tensor2d"
  // WARNING: all off-diagonal are discarded
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i )
      out[i] = m_data[ i*m_nd + i ];

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // norm
  // ----

  X norm() const { X C = static_cast<X>(0); for ( auto &i : m_data ) C += std::abs(i); return C; }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( auto &i : m_data ) i = static_cast<X>(0); };
  void ones        (     ) { for ( auto &i : m_data ) i = static_cast<X>(1); };
  void setConstant ( X D ) { for ( auto &i : m_data ) i = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2 <X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X           inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X> inline dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2 <X> inline T     (                    ) const; // transpose       : B_ij   = A_ji
  X           inline trace (                    ) const; // trace           : A_ii
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2 <X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i          )       { return m_data[i];        };
  const X& operator[](size_t i          ) const { return m_data[i];        };
  X&       operator()(size_t i, size_t j)       { return m_data[i*m_nd+j]; };
  const X& operator()(size_t i, size_t j) const { return m_data[i*m_nd+j]; };

  // arithmetic operators: tensor2 ?= tensor2
  // ----------------------------------------

  tensor2<X>& operator*= (const tensor2<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] *= B[i];

    return *this;
  };

  tensor2<X>& operator/= (const tensor2<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] /= B[i];

    return *this;
  };

  tensor2<X>& operator+= (const tensor2<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] += B[i];

    return *this;
  };

  tensor2<X>& operator-= (const tensor2<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor2 ?= tensor2s
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - extract value
        X b = B[ i*m_nd - (i-1)*i/2 + j - i ];
        // - store symmetrically
                      m_data[i*m_nd+j] *= b;
        if ( i != j ) m_data[j*m_nd+i] *= b;
      }
    }

    return *this;
  };

  tensor2<X>& operator/= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - extract value
        X b = B[ i*m_nd - (i-1)*i/2 + j - i ];
        // - store symmetrically
                      m_data[i*m_nd+j] /= b;
        if ( i != j ) m_data[j*m_nd+i] /= b;
      }
    }

    return *this;
  };

  tensor2<X>& operator+= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - extract value
        X b = B[ i*m_nd - (i-1)*i/2 + j - i ];
        // - store symmetrically
                      m_data[i*m_nd+j] += b;
        if ( i != j ) m_data[j*m_nd+i] += b;
      }
    }

    return *this;
  };

  tensor2<X>& operator-= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - extract value
        X b = B[ i*m_nd - (i-1)*i/2 + j - i ];
        // - store symmetrically
                      m_data[i*m_nd+j] -= b;
        if ( i != j ) m_data[j*m_nd+i] -= b;
      }
    }

    return *this;
  };

  // arithmetic operators: tensor2 ?= tensor2d
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=0; j<m_nd; ++j ) {
        if ( i == j ) m_data[i*m_nd+i] *= B[i];
        else          m_data[i*m_nd+j]  = static_cast<X>(0);
      }
    }

    return *this;
  };

  tensor2<X>& operator+= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i )
      m_data[i*m_nd+i] += B[i];

    return *this;
  };

  tensor2<X>& operator-= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i )
      m_data[i*m_nd+i] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor2 ?= scalar
  // ---------------------------------------

  tensor2<X>& operator*= (const X &B)
  {
    for ( auto &i: m_data )
      i *= B;

    return *this;
  };

  tensor2<X>& operator/= (const X &B)
  {
    for ( auto &i: m_data )
      i /= B;

    return *this;
  };

  tensor2<X>& operator+= (const X &B)
  {
    for ( auto &i: m_data )
      i += B;

    return *this;
  };

  tensor2<X>& operator-= (const X &B)
  {
    for ( auto &i: m_data )
      i -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0;  i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2s<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( m_data[i*m_nd+j] != B(i,j) )
          return false;

    return true;
  };

  bool operator== ( const tensor2d<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( m_data[i*m_nd+j] != B(i,j) )
          return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool issymmetric()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i+1 ; j < m_nd ; ++j )
        if ( m_data[ i*m_nd + j ] != m_data[ j*m_nd + i ] )
          return false;

    return true;
  };

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( i != j )
          if ( m_data[ i*m_nd + j ] )
            return false;

    return true;
  };

}; // class tensor2

// arithmetic operators: tensor2 = tensor2 ? tensor2
// -------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: tensor2 = tensor2 ? tensor2s
// --------------------------------------------------

template <class X> tensor2 <X> operator* (const tensor2 <X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] * b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] * b;
    }
  }

  return C;
}

template <class X> tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] / b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] / b;
    }
  }

  return C;
}

template <class X> tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] + b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] + b;
    }
  }

  return C;
}

template <class X> tensor2 <X> operator- (const tensor2 <X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] - b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] - b;
    }
  }

  return C;
}


// arithmetic operators: tensor2 = tensor2 ? tensor2d
// --------------------------------------------------

template <class X> tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i*nd + j ] + B[ i ];
      else          C[ i*nd + j ] = A[ i*nd + j ];
    }
  }

  return C;
}

template <class X> tensor2 <X> operator- (const tensor2 <X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i*nd + j ] - B[ i ];
      else          C[ i*nd + j ] = A[ i*nd + j ];
    }
  }

  return C;
}

// arithmetic operators: tensor2 = tensor2 ? scalar
// ------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: tensor2 = tensor2s ? tensor2
// --------------------------------------------------

template <class X> tensor2 <X> operator* (const tensor2s<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a * B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a * B[ j*nd + i ];
    }
  }

  return C;
}

template <class X> tensor2 <X> operator/ (const tensor2s<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a / B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a / B[ j*nd + i ];
    }
  }

  return C;
}

template <class X> tensor2 <X> operator+ (const tensor2s<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a + B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a + B[ j*nd + i ];
    }
  }

  return C;
}

template <class X> tensor2 <X> operator- (const tensor2s<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a - B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a - B[ j*nd + i ];
    }
  }

  return C;
}

// arithmetic operators: tensor2 = tensor2d ? tensor2
// --------------------------------------------------

template <class X> tensor2 <X> operator+ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i ] + B[ i*nd + j ];
      else          C[ i*nd + j ] =          B[ i*nd + j ];
    }
  }

  return C;
}

template <class X> tensor2 <X> operator- (const tensor2d<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2 <X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i ] - B[ i*nd + j ];
      else          C[ i*nd + j ] =        - B[ i*nd + j ];
    }
  }

  return C;
}

// arithmetic operators: tensor2 = scalar ? tensor2
// ------------------------------------------------

template <class X> tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A * B[i];

  return C;
}

template <class X> tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A / B[i];

  return C;
}

template <class X> tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A + B[i];

  return C;
}

template <class X> tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// cppmat::tensor2s (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor2s
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd;   // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2s(){};

  // explicit constructor: set correct size (WARNING: data not initialized)
  tensor2s(size_t nd ) { resize(nd); };

  // explicit constructor: set to constant "D"
  tensor2s(size_t nd, X D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i]=D; };

  // explicit constructor: from full matrix
  // WARNING: the input is symmetrized: "this(i,j) = ( in(i,j) + in(j,i) ) / 2."
  tensor2s(size_t nd, const X *D)
  {
    // check for symmetry (code eliminated if "NDEBUG" is defined at the beginning of the code)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < nd ; ++i )
        for ( size_t j = i+1 ; j < nd ; ++j )
          assert( D[ i*nd + j ] == D[ j*nd + i ] );
    #endif

    resize(nd);

    for ( size_t i = 0 ; i < nd ; ++i )
      for ( size_t j = i ; j < nd ; ++j )
        m_data[ i*nd - (i-1)*i/2 + j - i ] = D[ i*nd + j ];
  };

  // change number of dimensions (WARNING: data not initialized)
  void resize(size_t nd) { m_nd = nd; m_data.resize((nd+1)*nd/2); };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor2s", but of "cppmat::tensor2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,m_nd,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2s" -> "tensor2s" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2s<U> () const
  {
    tensor2s<U> out(m_nd);

    for ( size_t i=0; i<size(); ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy "const tensor2s" -> "tensor2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - get item
        U b = static_cast<U>( m_data[ i*m_nd - (i-1)*i/2 + j - i ] );
        // - store item, and symmetric copy
        out[ i*m_nd + j ] = b;
        out[ j*m_nd + i ] = b;
      }
    }

    return out;
  }

  // convert "tensor2s -> tensor2d"
  // WARNING: all off-diagonal are discarded
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i )
      out[i] = m_data[ i*m_nd - (i-1)*i/2 ];

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // norm
  // ----

  X norm() const { X C = static_cast<X>(0); for ( auto &i : m_data ) C += std::abs(i); return C; }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( auto &i : m_data ) i = static_cast<X>(0); };
  void ones        (     ) { for ( auto &i : m_data ) i = static_cast<X>(1); };
  void setConstant ( X D ) { for ( auto &i : m_data ) i = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2 <X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X           inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X> inline dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2s<X> inline T     (                    ) const; // transpose       : B_ij   = A_ji
  X           inline trace (                    ) const; // trace           : A_ii
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2s<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if (i <= j) return m_data[ i*m_nd - (i-1)*i/2 + j - i ];
    else        return m_data[ j*m_nd - (j-1)*j/2 + i - j ];
  }

  const X& operator()(size_t i, size_t j) const
  {
    if (i <= j) return m_data[ i*m_nd - (i-1)*i/2 + j - i ];
    else        return m_data[ j*m_nd - (j-1)*j/2 + i - j ];
  }

  // arithmetic operators: tensor2s ?= tensor2s
  // ------------------------------------------

  tensor2s<X>& operator*= (const tensor2s<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] *= B[i];

    return *this;
  };

  tensor2s<X>& operator/= (const tensor2s<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] /= B[i];

    return *this;
  };

  tensor2s<X>& operator+= (const tensor2s<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] += B[i];

    return *this;
  };

  tensor2s<X>& operator-= (const tensor2s<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor2s ?= tensor2d
  // ------------------------------------------

  tensor2s<X>& operator*= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        if ( i == j ) m_data[ i*m_nd - (i-1)*i/2         ] *= B[i];
        else          m_data[ i*m_nd - (i-1)*i/2 + j - i ]  = static_cast<X>(0);
      }
    }

    return *this;
  };

  tensor2s<X>& operator+= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i )
      m_data[ i*m_nd - (i-1)*i/2 ] += B[i];

    return *this;
  };

  tensor2s<X>& operator-= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i )
      m_data[ i*m_nd - (i-1)*i/2 ] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor2s ?= scalar
  // ----------------------------------------

  tensor2s<X>& operator*= (const X &B)
  {
    for ( auto &i: m_data )
      i *= B;

    return *this;
  };

  tensor2s<X>& operator/= (const X &B)
  {
    for ( auto &i: m_data )
      i /= B;

    return *this;
  };

  tensor2s<X>& operator+= (const X &B)
  {
    for ( auto &i: m_data )
      i += B;

    return *this;
  };

  tensor2s<X>& operator-= (const X &B)
  {
    for ( auto &i: m_data )
      i -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2s<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0; i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        if ( m_data[ i*m_nd - (i-1)*i/2 + j - i ] != B(i,j) ) return false;
        if ( m_data[ i*m_nd - (i-1)*i/2 + j - i ] != B(j,i) ) return false;
      }
    }

    return true;
  };

  bool operator== ( const tensor2d<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i ; j < m_nd ; ++j )
        if ( m_data[ i*m_nd - (i-1)*i/2 + j - i ] != B(i,j) ) return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i+1 ; j < m_nd ; ++j )
        if ( m_data[ i*m_nd - (i-1)*i/2 + j - i ] )
          return false;

    return true;
  };

}; // class tensor2s

// arithmetic operators: tensor2s = tensor2s ? tensor2s
// ----------------------------------------------------

template <class X> tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: tensor2s = tensor2s ? tensor2d
// ----------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i*nd - (i-1)*i/2         ] + B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i*nd - (i-1)*i/2         ] - B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// arithmetic operators: tensor2s = tensor2s ? scalar
// --------------------------------------------------

template <class X> tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X> tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: tensor2s = tensor2d ? scalar
// --------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] + B;
      else          C[ i*nd - (i-1)*i/2 + j - i ] =          B;
    }
  }

   return C;
}

template <class X> tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] - B;
      else          C[ i*nd - (i-1)*i/2 + j - i ] =        - B;
    }
  }

   return C;
}

// arithmetic operators: tensor2s = tensor2d ? tensor2s
// ----------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] + B[ i*nd - (i-1)*i/2         ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] =          B[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] - B[ i*nd - (i-1)*i/2         ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] =        - B[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// arithmetic operators: tensor2s = scalar ? tensor2s
// --------------------------------------------------

template <class X> tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A * B[i];

  return C;
}

template <class X> tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A / B[i];

  return C;
}

template <class X> tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A + B[i];

  return C;
}

template <class X> tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A - B[i];

  return C;
}

// arithmetic operators: tensor2s = scalar ? tensor2d
// --------------------------------------------------

template <class X> tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  size_t      nd = B.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A + B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A;
    }
  }

   return C;
}



template <class X> tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  size_t      nd = B.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A - B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A;
    }
  }

   return C;
}

// =================================================================================================
// cppmat::tensor2d (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor2d
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd;   // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2d(){};

  // explicit constructor: set correct size (WARNING: data not initialized)
  tensor2d(size_t nd ) { resize(nd); };

  // explicit constructor: set to constant "D"
  tensor2d(size_t nd, X D) { resize(nd); for ( size_t i=0; i<nd; ++i ) m_data[i]=D; };

  // explicit constructor: from full matrix
  // WARNING: all off-diagonal are discarded
  tensor2d(size_t nd, const X *D)
  {
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < nd ; ++i )
        for ( size_t j = 0 ; j < nd ; ++j )
          if ( i != j )
            assert( ! D[ i*nd + j ] );
    #endif

    resize(nd);

    for ( size_t i=0; i<nd; ++i )
      m_data[ i ] = D[ i*nd + i ];
  };

  // change number of dimensions (WARNING: data not initialized)
  // NB: a trick is used to transmit the off-diagonal zeros, therefore one zero is stored
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd+1); m_data[nd] = static_cast<X>(0); };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor2d", but of "cppmat::tensor2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,m_nd,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2d" -> "tensor2d" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2d<U> () const
  {
    tensor2d<U> out(m_nd);

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy "const tensor2d" -> "tensor2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd,static_cast<U>(0));

    for ( size_t i=0; i<m_nd; ++i )
      out[ i*m_nd + i ] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy "const tensor2d" -> "tensor2s" ( + change of type )
  // WARNING: all off-diagonal are discarded
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2s<U> () const
  {
    tensor2s<U> out(m_nd,static_cast<U>(0));

    for ( size_t i=0; i<m_nd; ++i )
      out[ i*m_nd - (i-1)*i/2 ] = static_cast<U>( m_data[i] );

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // norm
  // ----

  X norm() const { X C = static_cast<X>(0); for ( auto &i : m_data ) C += std::abs(i); return C; }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( size_t i=0; i<m_nd; ++i ) m_data[i] = static_cast<X>(0); };
  void ones        (     ) { for ( size_t i=0; i<m_nd; ++i ) m_data[i] = static_cast<X>(1); };
  void setConstant ( X D ) { for ( size_t i=0; i<m_nd; ++i ) m_data[i] = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2 <X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2d<X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X           inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X           inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X> inline dyadic(const tensor2 <X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2s<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor4 <X> inline dyadic(const tensor2d<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor2d<X> inline T     (                    ) const; // transpose      : B_ij   = A_ji
  X           inline trace (                    ) const; // trace          : A_ii
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2d<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if (i == j) return m_data[i];
    else        return m_data[m_nd];
  }

  const X& operator()(size_t i, size_t j) const
  {
    if (i == j) return m_data[i];
    else        return m_data[m_nd];
  }

  // arithmetic operators: tensor2d ?= tensor2d
  // ------------------------------------------

  tensor2d<X>& operator*= (const tensor2d<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] *= B[i];

    return *this;
  };

  tensor2d<X>& operator+= (const tensor2d<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] += B[i];

    return *this;
  };

  tensor2d<X>& operator-= (const tensor2d<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: tensor2d ?= tensor2
  // -----------------------------------------

  tensor2d<X>& operator*= (const tensor2 <X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] *= B[ i*m_nd+i ];

    return *this;
  };

  tensor2d<X>& operator/= (const tensor2 <X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] /= B[ i*m_nd+i ];

    return *this;
  };

  // arithmetic operators: tensor2d ?= tensor2s
  // ------------------------------------------

  tensor2d<X>& operator*= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] *= B[ i*m_nd - (i-1)*i/2 ];

    return *this;
  };

  tensor2d<X>& operator/= (const tensor2s<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] /= B[ i*m_nd - (i-1)*i/2 ];

    return *this;
  };

  // arithmetic operators: tensor2d ?= scalar
  // ----------------------------------------

  tensor2d<X>& operator*= (const X &B)
  {
    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] *= B;

    return *this;
  };

  tensor2d<X>& operator/= (const X &B)
  {
    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] /= B;

    return *this;
  };

  tensor2d<X>& operator+= (const X &B)
  {
    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] += B;

    return *this;
  };

  tensor2d<X>& operator-= (const X &B)
  {
    for ( size_t i=0; i<ndim(); ++i )
      m_data[i] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2d<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0;  i < m_nd ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = 0 ; j < m_nd ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

  bool operator== ( const tensor2s<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

}; // class tensor2d

// arithmetic operators: tensor2d = tensor2d ? tensor2d
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: tensor2d = tensor2d ? tensor2
// ---------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] * B[ i*nd + i ];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] / B[ i*nd + i ];

  return C;
}

// arithmetic operators: tensor2d = tensor2d ? tensor2s
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] * B[ i*nd - (i-1)*i/2 ];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] / B[ i*nd - (i-1)*i/2 ];

  return C;
}


// arithmetic operators: tensor2d = tensor2d ? scalar
// --------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A[i] / B;

  return C;
}

// arithmetic operators: tensor2d = tensor2 ? tensor2d
// ---------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[ i*nd + i ] * B[i];

  return C;
}

// arithmetic operators: tensor2d = tensor2s ? tensor2d
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t      nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[ i*nd - (i-1)*i/2 ] * B[i];

  return C;
}


// arithmetic operators: tensor2d = scalar ? tensor2d
// --------------------------------------------------

template <class X> tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  tensor2d<X> C(B.ndim());

  for ( size_t i=0; i<C.ndim(); ++i )
    C[i] = A * B[i];

  return C;
}

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X> class vector
{
private:

  std::vector<X> m_data; // data array
  size_t         m_nd;   // number of dimensions

public:

  // constructors
  // ------------

  // implicit constructor
  vector(){};

  // explicit constructor: set correct size (WARNING: data not initialized)
  vector(size_t nd ) { resize(nd); };

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  vector(size_t nd,       X  D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i] = D;    };
  vector(size_t nd, const X *D) { resize(nd); for ( size_t i=0; i<size(); ++i ) m_data[i] = D[i]; };

  // change number of dimensions (WARNING: data not initialized)
  void resize(size_t nd) { m_nd = nd; m_data.resize(nd); };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(1,m_nd,bytes); };

  // copy constructor
  // ----------------

  // copy "vector" -> "vector" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator vector<U> () const
  {
    vector<U> out(m_nd);

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return m_data.data (); };
  auto     begin()       { return m_data.begin(); };
  auto     end  ()       { return m_data.end  (); };

  // dimensions
  // ----------

  size_t size() const { return m_data.size(); };
  size_t ndim() const { return m_nd;          };

  // norm
  // ----

  X norm() const { X C = static_cast<X>(0); for ( auto &i : m_data ) C += std::abs(i); return C; }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( auto &i : m_data ) i = static_cast<X>(0); };
  void ones        (     ) { for ( auto &i : m_data ) i = static_cast<X>(1); };
  void setConstant ( X D ) { for ( auto &i : m_data ) i = D;                 };

  // tensor products / operations
  // ----------------------------

  X          inline dot   (const vector  <X> &B) const; // dot    product: C   = A_i*B_i
  vector <X> inline dot   (const tensor2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2<X> inline dyadic(const vector  <X> &B) const; // dyadic product: C_ij = A_i*B_j
  vector <X> inline cross (const vector  <X> &B) const; // cross product (only in 3D)

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; };
  const X& operator[](size_t i) const { return m_data[i]; };
  X&       operator()(size_t i)       { return m_data[i]; };
  const X& operator()(size_t i) const { return m_data[i]; };

  // arithmetic operators: vector ?= vector
  // --------------------------------------

  vector<X>& operator*= (const vector<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] *= B[i];

    return *this;
  };

  vector<X>& operator/= (const vector<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] /= B[i];

    return *this;
  };

  vector<X>& operator+= (const vector<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] += B[i];

    return *this;
  };

  vector<X>& operator-= (const vector<X> &B)
  {
    assert( size() == B.size() );
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<size(); ++i )
      m_data[i] -= B[i];

    return *this;
  };

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<X>& operator*= (const X &B)
  {
    for ( auto &i: m_data )
      i *= B;

    return *this;
  };

  vector<X>& operator/= (const X &B)
  {
    for ( auto &i: m_data )
      i /= B;

    return *this;
  };

  vector<X>& operator+= (const X &B)
  {
    for ( auto &i: m_data )
      i += B;

    return *this;
  };

  vector<X>& operator-= (const X &B)
  {
    for ( auto &i: m_data )
      i -= B;

    return *this;
  };

  // equality operator
  // -----------------

  bool operator== ( const vector<X> &B )
  {
    for ( size_t i = 0;  i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class vector

// arithmetic operators: vector = vector ? vector
  // --------------------------------------------

template <class X> vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X> vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: vector = vector ? scalar
  // --------------------------------------------

template <class X> vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X> vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X> vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X> vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
  // --------------------------------------------

template <class X> vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A * B[i];

  return C;
}

template <class X> vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A / B[i];

  return C;
}

template <class X> vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A + B[i];

  return C;
}

template <class X> vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i=0; i<C.size(); ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// print to screen
// =================================================================================================

template<class X> void inline tensor4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(m_nd).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l] );

}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor4<X>& src)
{
  for ( size_t i = 0 ; i < src.ndim() ; ++i )
    for ( size_t j = 0 ; j < src.ndim() ; ++j )
      for ( size_t k = 0 ; k < src.ndim() ; ++k )
        for ( size_t l = 0 ; l < src.ndim() ; ++l )
          out << "(" << i << "," << j << "," << k << "," << l << ") " << src(i,j,k,l) << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2<X>::printf(std::string fmt) const
{
  size_t nd = m_nd;

  for ( size_t h=0; h<nd; ++h ) {
    for ( size_t i=0; i<nd-1; ++i )
      std::printf((fmt+",").c_str(),m_data[h*nd+i]);
    std::printf((fmt+";\n").c_str(),m_data[h*nd+(nd-1)]);
  }
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2<X>& src)
{
  for ( size_t i=0; i<src.ndim(); ++i ) {
    for ( size_t j=0; j<src.ndim()-1; ++j )
      out << src(i,j) << ", ";
    out << src(i,src.ndim()-1) << "; " << std::endl;
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2s<X>::printf(std::string fmt) const
{
  size_t nd = m_nd;
  size_t i,j;

  for ( i=0; i<nd; ++i ) {
    for ( j=0; j<nd-1; ++j ) {
      if (i <= j) std::printf((fmt+",").c_str(),m_data[ i * nd - (i - 1) * i / 2 + j - i ]);
      else        std::printf((fmt+",").c_str(),m_data[ j * nd - (j - 1) * j / 2 + i - j ]);
    }
    j = nd-1;
    if (i <= j) std::printf((fmt+";\n").c_str(),m_data[ i * nd - (i - 1) * i / 2 + j - i ]);
    else        std::printf((fmt+";\n").c_str(),m_data[ j * nd - (j - 1) * j / 2 + i - j ]);
  }
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2s<X>& src)
{
  for ( size_t i=0; i<src.ndim(); ++i ) {
    for ( size_t j=0; j<src.ndim()-1; ++j )
      out << src(i,j) << ", ";
    out << src(i,src.ndim()-1) << "; " << std::endl;
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2d<X>::printf(std::string fmt) const
{
  size_t nd = m_nd;
  size_t i,j;

  for ( i=0; i<nd; ++i ) {
    for ( j=0; j<nd-1; ++j ) {
      if (i == j) std::printf((fmt+",").c_str(),m_data[i   ]);
      else        std::printf((fmt+",").c_str(),m_data[m_nd]);
    }
    j = nd-1;
    if (i == j) std::printf((fmt+";\n").c_str(),m_data[i   ]);
    else        std::printf((fmt+";\n").c_str(),m_data[m_nd]);
  }
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2d<X>& src)
{
  for ( size_t i=0; i<src.ndim(); ++i ) {
    for ( size_t j=0; j<src.ndim()-1; ++j )
      out << src(i,j) << ", ";
    out << src(i,src.ndim()-1) << "; " << std::endl;
  }
  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline vector<X>::printf(std::string fmt) const
{
  size_t nd = m_nd;

  for ( size_t i=0; i<nd-1; ++i )
    std::printf((fmt+",").c_str(),m_data[i]);
  std::printf((fmt+"\n").c_str(),m_data[nd-1]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  for ( size_t i=0; i<src.ndim()-1; ++i )
    out << src(i) << ", ";
  out << src(src.ndim()-1) << std::endl;
  return out;
}

// =================================================================================================
// identity tensors
// =================================================================================================

tensor4<double> inline identity4(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==l and j==k )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4rt(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==k and j==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4II(size_t nd)
{
  tensor4<double> I(nd,0.0);

  for ( size_t i=0; i<nd; ++i )
    for ( size_t j=0; j<nd; ++j )
      for ( size_t k=0; k<nd; ++k )
        for ( size_t l=0; l<nd; ++l )
          if ( i==j and k==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4s(size_t nd)
{ return (identity4(nd)+identity4rt(nd))/2.; }

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4d(size_t nd)
{ return identity4s(nd)-identity4II(nd)/static_cast<double>(nd); }

// -------------------------------------------------------------------------------------------------

tensor2d<double> inline identity2(size_t nd)
{
  tensor2d<double> I(nd,1.);

  return I;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X> tensor4<X> inline tensor4<X>::ddot(const tensor4<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          for ( size_t m=0; m<m_nd; ++m )
            for ( size_t n=0; n<m_nd; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l)*B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,j) += (*this)(i,j,k,k)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += (*this)(i,i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += ( m_data[i*m_nd-(i-1)*i/2] * B[i*m_nd-(i-1)*i/2] );

  for ( size_t i = 0 ; i<m_nd ; ++i )
    for ( size_t j = i+1 ; j<m_nd ; ++j )
      C += ( static_cast<X>(2) * m_data[i*m_nd-(i-1)*i/2+j-i] * B[i*m_nd-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += m_data[i*m_nd-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      for ( size_t l=0; l<m_nd; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += m_data[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += m_data[i]*B[i*m_nd-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += m_data[i]*B[i];

  return C;
}


// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2s<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::dot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2d<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    C(i,i) += (*this)(i,i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2d<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    C(i) += (*this)(i,i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    C(i) += (*this)(i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline vector<X>::dot(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += (*this)(i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      for ( size_t l=0; l<m_nd; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      for ( size_t l=0; l<m_nd; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t k=0; k<m_nd; ++k )
      C(i,i,k,k) += (*this)(i,i)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline vector<X>::dyadic(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C(m_nd,static_cast<X>(0));

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(i,j) += (*this)(i)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::cross(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  if ( m_nd != 3 )
    throw std::runtime_error("'cross' only implemented in 3D");

  vector<X> C(3);

  C[0] =                     m_data[1]*B[2]-B[1]*m_data[2] ;
  C[1] = static_cast<X>(-1)*(m_data[0]*B[2]-B[0]*m_data[2]);
  C[2] =                     m_data[0]*B[1]-B[0]*m_data[1] ;

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X> tensor4<X> inline tensor4<X>::T() const
{
  tensor4<X> C(m_nd);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::RT() const
{
  tensor4<X> C(m_nd);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::LT() const
{
  tensor4<X> C(m_nd);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      for ( size_t k=0; k<m_nd; ++k )
        for ( size_t l=0; l<m_nd; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::T() const
{
  tensor2<X> C(m_nd);

  for ( size_t i=0; i<m_nd; ++i )
    for ( size_t j=0; j<m_nd; ++j )
      C(j,i) = (*this)(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2s<X> inline tensor2s<X>::T() const
{
  tensor2s<X> C(m_nd);

  for ( size_t i=0; i<size(); ++i )
    C[i] = (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::T() const
{
  tensor2d<X> C(m_nd);

  for ( size_t i=0; i<size(); ++i )
    C[i] = (*this)[i];

  return C;
}

// =================================================================================================
// miscellaneous
// =================================================================================================

template<class X> X inline tensor2<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<m_nd; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<size(); ++i )
    C += (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::det() const
{
  if ( m_nd==2 )
   return m_data[0] * m_data[3] - m_data[1] * m_data[2];

  if ( m_nd==3 )
    return ( m_data[0] * m_data[4] * m_data[8] +
             m_data[1] * m_data[5] * m_data[6] +
             m_data[2] * m_data[3] * m_data[7] ) -
           ( m_data[2] * m_data[4] * m_data[6] +
             m_data[1] * m_data[3] * m_data[8] +
             m_data[0] * m_data[5] * m_data[7] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::det() const
{
  if ( m_nd==2 )
   return m_data[0] * m_data[2] - m_data[1] * m_data[1];

  if ( m_nd==3 )
    return (                     m_data[0] * m_data[3] * m_data[5] +
             static_cast<X>(2) * m_data[1] * m_data[2] * m_data[4] ) -
           (                     m_data[4] * m_data[4] * m_data[0] +
                                 m_data[2] * m_data[2] * m_data[3] +
                                 m_data[1] * m_data[1] * m_data[5] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::det() const
{
  X C = static_cast<X>(1);

  for ( size_t i=0; i<m_nd; ++i )
    C *= (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2<X> C(m_nd);

  if ( m_nd==2 ) {
    C[0] =                      m_data[3] / D;
    C[1] = static_cast<X>(-1) * m_data[1] / D;
    C[2] = static_cast<X>(-1) * m_data[2] / D;
    C[3] =                      m_data[0] / D;
    return C;
  }

  if ( m_nd==3 ) {
    C[0] = (m_data[4]*m_data[8]-m_data[5]*m_data[7]) / D;
    C[1] = (m_data[2]*m_data[7]-m_data[1]*m_data[8]) / D;
    C[2] = (m_data[1]*m_data[5]-m_data[2]*m_data[4]) / D;
    C[3] = (m_data[5]*m_data[6]-m_data[3]*m_data[8]) / D;
    C[4] = (m_data[0]*m_data[8]-m_data[2]*m_data[6]) / D;
    C[5] = (m_data[2]*m_data[3]-m_data[0]*m_data[5]) / D;
    C[6] = (m_data[3]*m_data[7]-m_data[4]*m_data[6]) / D;
    C[7] = (m_data[1]*m_data[6]-m_data[0]*m_data[7]) / D;
    C[8] = (m_data[0]*m_data[4]-m_data[1]*m_data[3]) / D;
    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2s<X> inline tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2s<X> C(m_nd);

  if ( m_nd==2 ) {
    C[0] =                      m_data[2] / D;
    C[1] = static_cast<X>(-1) * m_data[1] / D;
    C[2] =                      m_data[0] / D;
    return C;
  }

  if ( m_nd==3 ) {
    C[0] = (m_data[3]*m_data[5]-m_data[4]*m_data[4]) / D;
    C[1] = (m_data[2]*m_data[4]-m_data[1]*m_data[5]) / D;
    C[2] = (m_data[1]*m_data[4]-m_data[2]*m_data[3]) / D;
    C[3] = (m_data[0]*m_data[5]-m_data[2]*m_data[2]) / D;
    C[4] = (m_data[2]*m_data[1]-m_data[0]*m_data[4]) / D;
    C[5] = (m_data[0]*m_data[3]-m_data[1]*m_data[1]) / D;
    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C(m_nd);

  for ( size_t i = 0; i < m_nd ; ++i )
    C[i] = static_cast<X>(1) / m_data[i];

  return C;
}

// =================================================================================================
// create aliases to call class functions as functions, not members
// =================================================================================================

// products

// --

template<class X> tensor4 <X> inline ddot  (const tensor4 <X> &A, const tensor4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor4 <X> &A, const tensor2 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor4 <X> &A, const tensor2s<X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor4 <X> &A, const tensor2d<X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor2 <X> &A, const tensor4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor2s<X> &A, const tensor4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2 <X> inline ddot  (const tensor2d<X> &A, const tensor4 <X> &B)
{ return A.ddot(B); }

// --

template<class X>          X  inline ddot  (const tensor2 <X> &A, const tensor2 <X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2 <X> &A, const tensor2s<X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2 <X> &A, const tensor2d<X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2s<X> &A, const tensor2 <X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2s<X> &A, const tensor2s<X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2s<X> &A, const tensor2d<X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2d<X> &A, const tensor2 <X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2d<X> &A, const tensor2s<X> &B)
{ return A.ddot(B); }

template<class X>          X  inline ddot  (const tensor2d<X> &A, const tensor2d<X> &B)
{ return A.ddot(B); }

// --

template<class X> tensor2 <X> inline dot   (const tensor2 <X> &A, const tensor2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2 <X> &A, const tensor2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2 <X> &A, const tensor2d<X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2s<X> &A, const tensor2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2s<X> &A, const tensor2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2s<X> &A, const tensor2d<X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2d<X> &A, const tensor2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2 <X> inline dot   (const tensor2d<X> &A, const tensor2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2d<X> inline dot   (const tensor2d<X> &A, const tensor2d<X> &B)
{ return A.dot(B); }

// --

template<class X> vector  <X> inline dot   (const tensor2 <X> &A, const vector  <X> &B)
{ return A.dot(B); }

template<class X> vector  <X> inline dot   (const tensor2s<X> &A, const vector  <X> &B)
{ return A.dot(B); }

template<class X> vector  <X> inline dot   (const tensor2d<X> &A, const vector  <X> &B)
{ return A.dot(B); }

template<class X> vector  <X> inline dot   (const vector  <X> &A, const tensor2 <X> &B)
{ return A.dot(B); }

template<class X> vector  <X> inline dot   (const vector  <X> &A, const tensor2s<X> &B)
{ return A.dot(B); }

template<class X> vector  <X> inline dot   (const vector  <X> &A, const tensor2d<X> &B)
{ return A.dot(B); }

template<class X>          X  inline dot   (const vector  <X> &A, const vector  <X> &B)
{ return A.dot(B); }

// --

template<class X> tensor4 <X> inline dyadic(const tensor2 <X> &A, const tensor2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2 <X> &A, const tensor2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2 <X> &A, const tensor2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2s<X> &A, const tensor2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2s<X> &A, const tensor2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2s<X> &A, const tensor2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2d<X> &A, const tensor2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2d<X> &A, const tensor2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor4 <X> inline dyadic(const tensor2d<X> &A, const tensor2d<X> &B)
{ return A.dyadic(B); }

// --

template<class X> tensor2 <X> inline dyadic(const vector  <X> &A, const vector  <X> &B)
{ return A.dyadic(B); }

// --

template<class X> vector  <X> inline cross (const vector  <X> &A, const vector  <X> &B)
{ return A.cross (B); }

// operations
template<class X> tensor2 <X> inline transpose (const tensor2 <X> &A) { return A.T    (); }
template<class X> tensor2s<X> inline transpose (const tensor2s<X> &A) { return A.T    (); }
template<class X> tensor2d<X> inline transpose (const tensor2d<X> &A) { return A.T    (); }
template<class X> tensor4 <X> inline transpose (const tensor4 <X> &A) { return A.T    (); }
template<class X> tensor4 <X> inline transposeR(const tensor4 <X> &A) { return A.RT   (); }
template<class X> tensor4 <X> inline transposeL(const tensor4 <X> &A) { return A.LT   (); }
template<class X> tensor2 <X> inline inv       (const tensor2 <X> &A) { return A.inv  (); }
template<class X> tensor2s<X> inline inv       (const tensor2s<X> &A) { return A.inv  (); }
template<class X> tensor2d<X> inline inv       (const tensor2d<X> &A) { return A.inv  (); }
template<class X>          X  inline det       (const tensor2 <X> &A) { return A.det  (); }
template<class X>          X  inline det       (const tensor2s<X> &A) { return A.det  (); }
template<class X>          X  inline det       (const tensor2d<X> &A) { return A.det  (); }
template<class X>          X  inline trace     (const tensor2 <X> &A) { return A.trace(); }
template<class X>          X  inline trace     (const tensor2s<X> &A) { return A.trace(); }
template<class X>          X  inline trace     (const tensor2d<X> &A) { return A.trace(); }

// =================================================================================================

} // namespace tensor

#endif

