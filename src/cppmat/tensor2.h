/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef TENSOR2_H
#define TENSOR2_H

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#include "tensor.h"

namespace cppmat {

// =================================================================================================
// forward declaration
// =================================================================================================

template<class X> class tensor2_4;
template<class X> class tensor2_2;
template<class X> class tensor2_2s;
template<class X> class tensor2_2d;
template<class X> class vector2;

// =================================================================================================
// cppmat::tensor2_4
// =================================================================================================

template<class X> class tensor2_4
{
private:

  X m_data[16]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2_4(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor2_4(      X  D) { for ( size_t i=0; i<16; ++i ) m_data[i]=D;    };
  tensor2_4(const X *D) { for ( size_t i=0; i<16; ++i ) m_data[i]=D[i]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(4,2,bytes); };

  // copy constructor
  // ----------------

  // copy "tensor2_4" -> "tensor2_4" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_4<U> () const
  {
    tensor2_4<U> out;

    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
      out[i+3] = static_cast<U>( m_data[i+3] );
    }

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // pointer
  // -------

  const X* data () const { return &m_data[0]; };

  // dimensions
  // ----------

  size_t size() const { return 16; };
  size_t ndim() const { return 3;  };

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      C += std::abs(m_data[i  ]);
      C += std::abs(m_data[i+1]);
      C += std::abs(m_data[i+2]);
      C += std::abs(m_data[i+3]);
    }

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( size_t i=0; i<16; ++i ) m_data[i] = static_cast<X>(0); }
  void ones        (     ) { for ( size_t i=0; i<16; ++i ) m_data[i] = static_cast<X>(1); }
  void setConstant ( X D ) { for ( size_t i=0; i<16; ++i ) m_data[i] = D;                 }

  // tensor products / operations
  // ----------------------------

  tensor2_4 <X> inline ddot(const tensor2_4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2_2 <X> inline ddot(const tensor2_2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2_2 <X> inline ddot(const tensor2_2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2_2 <X> inline ddot(const tensor2_2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2_4 <X> inline T   (                      ) const; // transposition   : B_lkji = A_ijkl
  tensor2_4 <X> inline RT  (                      ) const; // transposition   : B_ijlk = A_ijkl
  tensor2_4 <X> inline LT  (                      ) const; // transposition   : B_jikl = A_ijkl

  // index operators
  // ---------------

  X& operator[](size_t i)
  { return m_data[i]; };

  const X& operator[](size_t i) const
  { return m_data[i]; };

  X& operator()(size_t i, size_t j, size_t k, size_t l)
  { return m_data[i*8+j*4+k*2+l]; };

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*8+j*4+k*2+l]; };

  // arithmetic operators: tensor2_4 ?= tensor2_4
  // --------------------------------------------

  tensor2_4<X>& operator*= (const tensor2_4<X> &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] *= B[i  ];
      m_data[i+1] *= B[i+1];
      m_data[i+2] *= B[i+2];
      m_data[i+3] *= B[i+3];
    }

    return *this;
  };

  tensor2_4<X>& operator/= (const tensor2_4<X> &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] /= B[i  ];
      m_data[i+1] /= B[i+1];
      m_data[i+2] /= B[i+2];
      m_data[i+3] /= B[i+3];
    }

    return *this;
  };

  tensor2_4<X>& operator+= (const tensor2_4<X> &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] += B[i  ];
      m_data[i+1] += B[i+1];
      m_data[i+2] += B[i+2];
      m_data[i+3] += B[i+3];
    }

    return *this;
  };

  tensor2_4<X>& operator-= (const tensor2_4<X> &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] -= B[i  ];
      m_data[i+1] -= B[i+1];
      m_data[i+2] -= B[i+2];
      m_data[i+3] -= B[i+3];
    }

    return *this;
  };

  // arithmetic operators: tensor2_4 ?= tensor2_4
  // --------------------------------------------

  tensor2_4<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] *= B;
      m_data[i+1] *= B;
      m_data[i+2] *= B;
      m_data[i+3] *= B;
    }

    return *this;
  };

  tensor2_4<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] /= B;
      m_data[i+1] /= B;
      m_data[i+2] /= B;
      m_data[i+3] /= B;
    }

    return *this;
  };

  tensor2_4<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] += B;
      m_data[i+1] += B;
      m_data[i+2] += B;
      m_data[i+3] += B;
    }

    return *this;
  };

  tensor2_4<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < 16 ; i += 4 ) {
      m_data[i  ] -= B;
      m_data[i+1] -= B;
      m_data[i+2] -= B;
      m_data[i+3] -= B;
    }

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2_4<X> &B )
  {
    for ( size_t i = 0 ; i<16 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class tensor2_4

// arithmetic operators: tensor2_4 = tensor2_4 ? tensor2_4
// -------------------------------------------------------

template <class X> tensor2_4<X> operator* (const tensor2_4<X> &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
    C[i+3] = A[i+3] * B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator/ (const tensor2_4<X> &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
    C[i+3] = A[i+3] / B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator+ (const tensor2_4<X> &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
    C[i+3] = A[i+3] + B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator- (const tensor2_4<X> &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
    C[i+3] = A[i+3] - B[i+3];
  }

  return C;
}

// arithmetic operators: tensor2_4 = tensor2_4 ? scalar
// ----------------------------------------------------

template <class X> tensor2_4<X> operator* (const tensor2_4<X> &A, const X &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
    C[i+3] = A[i+3] * B;
  }

  return C; }

template <class X> tensor2_4<X> operator/ (const tensor2_4<X> &A, const X &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
    C[i+3] = A[i+3] / B;
  }

  return C; }

template <class X> tensor2_4<X> operator+ (const tensor2_4<X> &A, const X &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
    C[i+3] = A[i+3] + B;
  }

  return C; }

template <class X> tensor2_4<X> operator- (const tensor2_4<X> &A, const X &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
    C[i+3] = A[i+3] - B;
  }

  return C; }

// arithmetic operators: tensor2_4 = scalar ? tensor2_4
// ----------------------------------------------------

template <class X> tensor2_4<X> operator* (const X &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
    C[i+3] = A * B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator/ (const X &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
    C[i+3] = A / B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator+ (const X &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
    C[i+3] = A + B[i+3];
  }

  return C;
}

template <class X> tensor2_4<X> operator- (const X &A, const tensor2_4<X> &B)
{
  tensor2_4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 ) {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
    C[i+3] = A - B[i+3];
  }

  return C;
}

// =================================================================================================
// cppmat::tensor2_2
// =================================================================================================

template<class X> class tensor2_2
{
private:

  X m_data[4]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2_2(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor2_2(      X  D) { for ( size_t i=0; i<4; ++i ) m_data[i]=D;    };
  tensor2_2(const X *D) { for ( size_t i=0; i<4; ++i ) m_data[i]=D[i]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,2,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2_2" -> "tensor2_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2<U> () const
  {
    tensor2_2<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );
    out[3] = static_cast<U>( m_data[3] );

    return out;
  }

  // convert "tensor2_2 -> tensor2_2s"
  // WARNING: the output is symmetrized: "out(i,j) = ( this(i,j) + this(j,i) ) / 2."
  tensor2_2s<X> astensor2s()
  {
    tensor2_2s<X> out;

    out[0] =   m_data[0];
    out[1] = ( m_data[1] + m_data[2] ) / static_cast<X>(2);
    out[2] =   m_data[3];

    return out;
  }

  // convert "tensor2_2 -> tensor2_2d"
  // WARNING: all off-diagonal are discarded
  tensor2_2d<X> astensor2d()
  {
    tensor2_2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[3];

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return &m_data[0]; };

  // dimensions
  // ----------

  size_t size() const { return 4; };
  size_t ndim() const { return 2; };

  // norm
  // ----

  X norm() const
  {
    X C;

    C  = std::abs(m_data[0]);
    C += std::abs(m_data[1]);
    C += std::abs(m_data[2]);
    C += std::abs(m_data[3]);

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( size_t i=0; i<4; ++i ) m_data[i] = static_cast<X>(0); };
  void ones        (     ) { for ( size_t i=0; i<4; ++i ) m_data[i] = static_cast<X>(1); };
  void setConstant ( X D ) { for ( size_t i=0; i<4; ++i ) m_data[i] = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2_2 <X> inline dot   (const tensor2_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2 <X> inline dot   (const tensor2_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2 <X> inline dot   (const tensor2_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector2   <X> inline dot   (const vector2   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2_2 <X> inline ddot  (const tensor2_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor2_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor2_4 <X> inline dyadic(const tensor2_2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_2 <X> inline T     (                      ) const; // transpose       : B_ij   = A_ji
  X             inline trace (                      ) const; // trace           : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor2_2 <X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i          )       { return m_data[i];     };
  const X& operator[](size_t i          ) const { return m_data[i];     };
  X&       operator()(size_t i, size_t j)       { return m_data[i*2+j]; };
  const X& operator()(size_t i, size_t j) const { return m_data[i*2+j]; };

  // arithmetic operators: tensor2_2 ?= tensor2_2
  // --------------------------------------------

  tensor2_2<X>& operator*= (const tensor2_2<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];
    m_data[3] *= B[3];

    return *this;
  };

  tensor2_2<X>& operator/= (const tensor2_2<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];
    m_data[3] /= B[3];

    return *this;
  };

  tensor2_2<X>& operator+= (const tensor2_2<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];
    m_data[3] += B[3];

    return *this;
  };

  tensor2_2<X>& operator-= (const tensor2_2<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];
    m_data[3] -= B[3];

    return *this;
  };

  // arithmetic operators: tensor2_2 ?= tensor2_2s
  // ---------------------------------------------

  tensor2_2<X>& operator*= (const tensor2_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1]; m_data[2] *= B[1];
    m_data[3] *= B[2];

    return *this;
  };

  tensor2_2<X>& operator/= (const tensor2_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1]; m_data[2] /= B[1];
    m_data[3] /= B[2];

    return *this;
  };

  tensor2_2<X>& operator+= (const tensor2_2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1]; m_data[2] += B[1];
    m_data[3] += B[2];

    return *this;
  };

  tensor2_2<X>& operator-= (const tensor2_2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1]; m_data[2] -= B[1];
    m_data[3] -= B[2];

    return *this;
  };

  // arithmetic operators: tensor2_2 ?= tensor2_2d
  // -----------------------------------------

  tensor2_2<X>& operator*= (const tensor2_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[3] *= B[1];
    m_data[1] = m_data[2] = static_cast<X>(0);

    return *this;
  };

  tensor2_2<X>& operator+= (const tensor2_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[3] += B[1];

    return *this;
  };

  tensor2_2<X>& operator-= (const tensor2_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[3] -= B[1];

    return *this;
  };

  // arithmetic operators: tensor2_2 ?= scalar
  // ---------------------------------------

  tensor2_2<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;
    m_data[3] *= B;

    return *this;
  };

  tensor2_2<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;
    m_data[3] /= B;

    return *this;
  };

  tensor2_2<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;
    m_data[3] += B;

    return *this;
  };

  tensor2_2<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;
    m_data[3] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2_2<X> &B )
  {
    for ( size_t i = 0;  i < 4 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2_2s<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = 0 ; j < 2 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  };

  bool operator== ( const tensor2_2d<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = 0 ; j < 2 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool issymmetric()
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = i+1 ; j < 2 ; ++j )
        if ( m_data[ i*3 + j ] != m_data[ j*3 + i ] )
          return false;

    return true;
  };

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = 0 ; j < 2 ; ++j )
        if ( i != j )
          if ( m_data[ i*3 + j ] )
            return false;

    return true;
  };

}; // class tensor2_2

// arithmetic operators: tensor2_2 = tensor2_2 ? tensor2_2
// -------------------------------------------------------

template <class X> tensor2_2<X> operator* (const tensor2_2<X> &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];

  return C;
}

template <class X> tensor2_2<X> operator/ (const tensor2_2<X> &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];

  return C;
}

template <class X> tensor2_2<X> operator+ (const tensor2_2<X> &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];

  return C;
}

template <class X> tensor2_2<X> operator- (const tensor2_2<X> &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];

  return C;
}

// arithmetic operators: tensor2_2 = tensor2_2 ? tensor2_2s
// --------------------------------------------------------

template <class X> tensor2_2 <X> operator* (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[2] * B[1];
  C[3] = A[3] * B[2];

  return C;
}

template <class X> tensor2_2 <X> operator/ (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[2] / B[1];
  C[3] = A[3] / B[2];

  return C;
}

template <class X> tensor2_2 <X> operator+ (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[2] + B[1];
  C[3] = A[3] + B[2];

  return C;
}

template <class X> tensor2_2 <X> operator- (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[2] - B[1];
  C[3] = A[3] - B[2];

  return C;
}


// arithmetic operators: tensor2_2 = tensor2_2 ? tensor2_2d
// --------------------------------------------------------

template <class X> tensor2_2 <X> operator+ (const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];

  return C;
}

template <class X> tensor2_2 <X> operator- (const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];

  return C;
}

// arithmetic operators: tensor2_2 = tensor2_2 ? scalar
// ----------------------------------------------------

template <class X> tensor2_2<X> operator* (const tensor2_2<X> &A, const X &B)
{
  tensor2_2<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;

  return C;
}

template <class X> tensor2_2<X> operator/ (const tensor2_2<X> &A, const X &B)
{
  tensor2_2<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;

  return C;
}

template <class X> tensor2_2<X> operator+ (const tensor2_2<X> &A, const X &B)
{
  tensor2_2<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;

  return C;
}

template <class X> tensor2_2<X> operator- (const tensor2_2<X> &A, const X &B)
{
  tensor2_2<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;

  return C;
}

// arithmetic operators: tensor2_2 = tensor2_2s ? tensor2_2
// --------------------------------------------------

template <class X> tensor2_2 <X> operator* (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[1] * B[2];
  C[3] = A[2] * B[3];

  return C;
}

template <class X> tensor2_2 <X> operator/ (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[1] / B[2];
  C[3] = A[2] / B[3];

  return C;
}

template <class X> tensor2_2 <X> operator+ (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[1] + B[2];
  C[3] = A[2] + B[3];

  return C;
}

template <class X> tensor2_2 <X> operator- (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[1] - B[2];
  C[3] = A[2] - B[3];

  return C;
}

// arithmetic operators: tensor2_2 = tensor2_2d ? tensor2_2
// --------------------------------------------------

template <class X> tensor2_2 <X> operator+ (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];

  return C;
}

template <class X> tensor2_2 <X> operator- (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];

  return C;
}

// arithmetic operators: tensor2_2 = scalar ? tensor2_2
// ------------------------------------------------

template <class X> tensor2_2<X> operator* (const X &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];

  return C;
}

template <class X> tensor2_2<X> operator/ (const X &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];

  return C;
}

template <class X> tensor2_2<X> operator+ (const X &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];

  return C;
}

template <class X> tensor2_2<X> operator- (const X &A, const tensor2_2<X> &B)
{
  tensor2_2<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];

  return C;
}

// =================================================================================================
// cppmat::tensor2_2s (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor2_2s
{
private:

  X m_data[3]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2_2s(){};

  // explicit constructor: set to constant "D"
  tensor2_2s(X D) { for ( size_t i=0; i<3; ++i ) m_data[i]=D; };

  // explicit constructor: from full matrix
  // WARNING: the input is symmetrized: "this(i,j) = ( in(i,j) + in(j,i) ) / 2."
  tensor2_2s(const X *D)
  {
    // check for symmetry (code eliminated if "NDEBUG" is defined at the beginning of the code)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 2 ; ++i )
        for ( size_t j = i+1 ; j < 2 ; ++j )
          assert( D[ i*2 + j ] == D[ j*2 + i ] );
    #endif

    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = i ; j < 2 ; ++j )
        m_data[ i*2 - (i-1)*i/2 + j - i ] = D[ i*2 + j ];
  };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor2_2s", but of "cppmat::tensor2_2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,2,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2_2s" -> "tensor2_2s" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2s<U> () const
  {
    tensor2_2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // copy "const tensor2_2s" -> "tensor2_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2<U> () const
  {
    tensor2_2<U> out;

    out[0]          = m_data[0];
    out[1] = out[2] = m_data[1];
    out[3]          = m_data[2];

    return out;
  }

  // convert "tensor2_2s -> tensor2_2d"
  // WARNING: all off-diagonal are discarded
  tensor2_2d<X> astensor2d()
  {
    tensor2_2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[2];

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return &m_data[0]; };

  // dimensions
  // ----------

  size_t size() const { return 3; };
  size_t ndim() const { return 2; };

  // norm
  // ----

  X norm() const
  {
    X C;

    C  = std::abs(m_data[0]);
    C += std::abs(m_data[1]);
    C += std::abs(m_data[2]);

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros()             { m_data[0]=m_data[1]=m_data[2]=static_cast<X>(0); };
  void ones ()             { m_data[0]=m_data[1]=m_data[2]=static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0]=m_data[1]=m_data[2]=D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2_2 <X> inline dot   (const tensor2_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2 <X> inline dot   (const tensor2_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2 <X> inline dot   (const tensor2_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector2   <X> inline dot   (const vector2   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2_2 <X> inline ddot  (const tensor2_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor2_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor2_4 <X> inline dyadic(const tensor2_2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2_2s<X> inline T     (                      ) const; // transpose       : B_ij   = A_ji
  X             inline trace (                      ) const; // trace           : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor2_2s<X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if ( i == 0 ) {
      if ( j == 0 ) return m_data[0];
      else          return m_data[1];
    }
    else {
      if ( j == 0 ) return m_data[1];
      else          return m_data[2];
    }
  }

  const X& operator()(size_t i, size_t j) const
  {
    if ( i == 0 ) {
      if ( j == 0 ) return m_data[0];
      else          return m_data[1];
    }
    else {
      if ( j == 0 ) return m_data[1];
      else          return m_data[2];
    }
  }

  // arithmetic operators: tensor2_2s ?= tensor2_2s
  // ------------------------------------------

  tensor2_2s<X>& operator*= (const tensor2_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];

    return *this;
  };

  tensor2_2s<X>& operator/= (const tensor2_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];

    return *this;
  };

  tensor2_2s<X>& operator+= (const tensor2_2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];

    return *this;
  };

  tensor2_2s<X>& operator-= (const tensor2_2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];

    return *this;
  };

  // arithmetic operators: tensor2_2s ?= tensor2_2d
  // ------------------------------------------

  tensor2_2s<X>& operator*= (const tensor2_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[2] *= B[1];
    m_data[1]  = static_cast<X>(0);

    return *this;
  };

  tensor2_2s<X>& operator+= (const tensor2_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[2] += B[1];

    return *this;
  };

  tensor2_2s<X>& operator-= (const tensor2_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[2] -= B[1];

    return *this;
  };

  // arithmetic operators: tensor2_2s ?= scalar
  // ----------------------------------------

  tensor2_2s<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;

    return *this;
  };

  tensor2_2s<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;

    return *this;
  };

  tensor2_2s<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;

    return *this;
  };

  tensor2_2s<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2_2s<X> &B )
  {
    for ( size_t i = 0; i < 2 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2_2<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i ) {
      for ( size_t j = i ; j < 2 ; ++j ) {
        if ( m_data[ i*2 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;
        if ( m_data[ i*2 - (i-1)*i/2 + j - i ] != B(j,i) ) return false;
      }
    }

    return true;
  };

  bool operator== ( const tensor2_2d<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = i ; j < 2 ; ++j )
        if ( m_data[ i*2 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 2 ; ++i )
      for ( size_t j = i+1 ; j < 2 ; ++j )
        if ( m_data[ i*2 - (i-1)*i/2 + j - i ] )
          return false;

    return true;
  };

}; // class tensor2_2s

// arithmetic operators: tensor2_2s = tensor2_2s ? tensor2_2s
// ----------------------------------------------------------

template <class X> tensor2_2s<X> operator* (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

template <class X> tensor2_2s<X> operator/ (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

template <class X> tensor2_2s<X> operator+ (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

template <class X> tensor2_2s<X> operator- (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// arithmetic operators: tensor2_2s = tensor2_2s ? tensor2_2d
// ----------------------------------------------------------

template <class X> tensor2_2s<X> operator+ (const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0]+B[0];
  C[1] = A[1];
  C[2] = A[2]+B[1];

  return C;
}

template <class X> tensor2_2s<X> operator- (const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0]-B[0];
  C[1] = A[1];
  C[2] = A[2]-B[1];

  return C;
}

// arithmetic operators: tensor2_2s = tensor2_2s ? scalar
// ------------------------------------------------------

template <class X> tensor2_2s<X> operator* (const tensor2_2s<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

template <class X> tensor2_2s<X> operator/ (const tensor2_2s<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

template <class X> tensor2_2s<X> operator+ (const tensor2_2s<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

template <class X> tensor2_2s<X> operator- (const tensor2_2s<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// arithmetic operators: tensor2_2s = tensor2_2d ? scalar
// ------------------------------------------------------

template <class X> tensor2_2s<X> operator+ (const tensor2_2d<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] = A[1] + B;

  return C;
}

template <class X> tensor2_2s<X> operator- (const tensor2_2d<X> &A, const X &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] = A[1] - B;

  return C;
}

// arithmetic operators: tensor2_2s = tensor2_2d ? tensor2_2s
// ----------------------------------------------------------

template <class X> tensor2_2s<X> operator+ (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] = A[1] + B[2];

  return C;
}

template <class X> tensor2_2s<X> operator- (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] = A[1] - B[2];

  return C;
}

// arithmetic operators: tensor2_2s = scalar ? tensor2_2s
// ------------------------------------------------------

template <class X> tensor2_2s<X> operator* (const X &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

template <class X> tensor2_2s<X> operator/ (const X &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

template <class X> tensor2_2s<X> operator+ (const X &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

template <class X> tensor2_2s<X> operator- (const X &A, const tensor2_2s<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// arithmetic operators: tensor2_2s = scalar ? tensor2_2d
// ------------------------------------------------------

template <class X> tensor2_2s<X> operator+ (const X &A, const tensor2_2d<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A + B[1];

  return C;
}

template <class X> tensor2_2s<X> operator- (const X &A, const tensor2_2d<X> &B)
{
  tensor2_2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A - B[1];

  return C;
}

// =================================================================================================
// cppmat::tensor2_2d (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor2_2d
{
private:

  X m_data[3]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor2_2d(){ m_data[2] = static_cast<X>(0); };

  // explicit constructor: set to constant "D"
  tensor2_2d(X D) { m_data[0] = m_data[1] = D; };

  // explicit constructor: from full matrix
  // WARNING: all off-diagonal are discarded
  tensor2_2d(const X *D)
  {
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 2 ; ++i )
        for ( size_t j = 0 ; j < 2 ; ++j )
          if ( i != j )
            assert( ! D[ i*2 + j ] );
    #endif

    for ( size_t i=0; i<2; ++i )
      m_data[ i ] = D[ i*2 + i ];
  };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor2_2d", but of "cppmat::tensor2_2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,2,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor2_2d" -> "tensor2_2d" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2d<U> () const
  {
    tensor2_2d<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );

    return out;
  }

  // copy "const tensor2_2d" -> "tensor2_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2<U> () const
  {
    tensor2_2<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[3] = static_cast<U>( m_data[1] );

    out[1] = out[2] = static_cast<U>(0);

    return out;
  }

  // copy "const tensor2_2d" -> "tensor2_2s" ( + change of type )
  // WARNING: all off-diagonal are discarded
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor2_2s<U> () const
  {
    tensor2_2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[2] = static_cast<U>( m_data[1] );

    out[1] = static_cast<U>(0);

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return &m_data[0]; };

  // dimensions
  // ----------

  size_t size() const { return 3; };
  size_t ndim() const { return 2; };

  // norm
  // ----

  X norm() const
  {
    X C;

    C  = std::abs(m_data[0]);
    C += std::abs(m_data[1]);

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { m_data[0] = m_data[1] = static_cast<X>(0); };
  void ones        (     ) { m_data[0] = m_data[1] = static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0] = m_data[1] = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor2_2 <X> inline dot   (const tensor2_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2 <X> inline dot   (const tensor2_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2_2d<X> inline dot   (const tensor2_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector2   <X> inline dot   (const vector2   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2_2 <X> inline ddot  (const tensor2_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor2_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor2_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor2_4 <X> inline dyadic(const tensor2_2 <X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2s<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor2_4 <X> inline dyadic(const tensor2_2d<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor2_2d<X> inline T     (                      ) const; // transpose      : B_ij   = A_ji
  X             inline trace (                      ) const; // trace          : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor2_2d<X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if (i == j) return m_data[i];
    else        return m_data[2];
  }

  const X& operator()(size_t i, size_t j) const
  {
    if (i == j) return m_data[i];
    else        return m_data[2];
  }

  // arithmetic operators: tensor2_2d ?= tensor2_2d
  // ----------------------------------------------

  tensor2_2d<X>& operator*= (const tensor2_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];

    return *this;
  };

  tensor2_2d<X>& operator+= (const tensor2_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];

    return *this;
  };

  tensor2_2d<X>& operator-= (const tensor2_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];

    return *this;
  };

  // arithmetic operators: tensor2_2d ?= tensor2_2
  // ---------------------------------------------

  tensor2_2d<X>& operator*= (const tensor2_2 <X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[3];

    return *this;
  };

  tensor2_2d<X>& operator/= (const tensor2_2 <X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[3];

    return *this;
  };

  // arithmetic operators: tensor2_2d ?= tensor2_2s
  // ----------------------------------------------

  tensor2_2d<X>& operator*= (const tensor2_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[2];

    return *this;
  };

  tensor2_2d<X>& operator/= (const tensor2_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[2];

    return *this;
  };

  // arithmetic operators: tensor2_2d ?= scalar
  // ------------------------------------------

  tensor2_2d<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;

    return *this;
  };

  tensor2_2d<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;

    return *this;
  };

  tensor2_2d<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;

    return *this;
  };

  tensor2_2d<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor2_2d<X> &B )
  {
    for ( size_t i = 0;  i < 2 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor2_2<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i ) {
      for ( size_t j = 0 ; j < 2 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

  bool operator== ( const tensor2_2s<X> &B )
  {
    for ( size_t i = 0 ; i < 2 ; ++i ) {
      for ( size_t j = i ; j < 2 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

}; // class tensor2_2d

// arithmetic operators: tensor2_2d = tensor2_2d ? tensor2_2d
// ----------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

template <class X> tensor2_2d<X> operator+ (const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

template <class X> tensor2_2d<X> operator- (const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// arithmetic operators: tensor2_2d = tensor2_2d ? tensor2_2
// ---------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];

  return C;
}

template <class X> tensor2_2d<X> operator/ (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];

  return C;
}

// arithmetic operators: tensor2_2d = tensor2_2d ? tensor2_2s
// ----------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[2];

  return C;
}

template <class X> tensor2_2d<X> operator/ (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[2];

  return C;
}

// arithmetic operators: tensor2_2d = tensor2_2d ? scalar
// ------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2d<X> &A, const X &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

template <class X> tensor2_2d<X> operator/ (const tensor2_2d<X> &A, const X &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// arithmetic operators: tensor2_2d = tensor2_2 ? tensor2_2d
// ---------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];

  return C;
}

// arithmetic operators: tensor2_2d = tensor2_2s ? tensor2_2d
// ----------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[2] * B[1];

  return C;
}


// arithmetic operators: tensor2_2d = scalar ? tensor2_2d
// ------------------------------------------------------

template <class X> tensor2_2d<X> operator* (const X &A, const tensor2_2d<X> &B)
{
  tensor2_2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// =================================================================================================
// cppmat::vector2
// =================================================================================================

template<class X> class vector2
{
private:

  X m_data[2]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  vector2(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  vector2(      X  D) { m_data[0] = m_data[1] = D; };
  vector2(const X *D) { m_data[0] = D[0]; m_data[1] = D[1]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(1,2,bytes); };

  // copy constructor
  // ----------------

  // copy "vector2" -> "vector2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator vector2<U> () const
  {
    vector2<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );

    return out;
  }

  // print to screen
  // ---------------

  // formatted print (code below); NB "operator<<" is defined below
  void printf(std::string fmt) const;

  // iterators / pointer
  // -------------------

  const X* data () const { return &m_data[0]; };

  // dimensions
  // ----------

  size_t size() const { return 2; };
  size_t ndim() const { return 2; };

  // norm
  // ----

  X norm() const
  {
    X C;

    C  = std::abs(m_data[0]);
    C += std::abs(m_data[1]);

    return C;
  }

  // length
  // ------

  X length() const
  {
    X C;

    C  = std::pow(m_data[0],2.);
    C += std::pow(m_data[1],2.);

    return std::pow(C,.5);
  }

  // normalize to unit length
  // ------------------------

  void setUnitLength()
  {
    X C = length();

    if ( C <= static_cast<X>(0) ) return;

    m_data[0] /= C;
    m_data[1] /= C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { m_data[0] = m_data[1] = static_cast<X>(0); };
  void ones        (     ) { m_data[0] = m_data[1] = static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0] = m_data[1] = D;                 };

  // tensor products / operations
  // ----------------------------

  X            inline dot   (const vector2   <X> &B) const; // dot    product: C   = A_i*B_i
  vector2  <X> inline dot   (const tensor2_2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector2  <X> inline dot   (const tensor2_2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector2  <X> inline dot   (const tensor2_2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2_2<X> inline dyadic(const vector2   <X> &B) const; // dyadic product: C_ij = A_i*B_j

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; };
  const X& operator[](size_t i) const { return m_data[i]; };
  X&       operator()(size_t i)       { return m_data[i]; };
  const X& operator()(size_t i) const { return m_data[i]; };

  // arithmetic operators: vector2 ?= vector2
  // --------------------------------------

  vector2<X>& operator*= (const vector2<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];

    return *this;
  };

  vector2<X>& operator/= (const vector2<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];

    return *this;
  };

  vector2<X>& operator+= (const vector2<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];

    return *this;
  };

  vector2<X>& operator-= (const vector2<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];

    return *this;
  };

  // arithmetic operators: vector2 ?= scalar
  // --------------------------------------

  vector2<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;

    return *this;
  };

  vector2<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;

    return *this;
  };

  vector2<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;

    return *this;
  };

  vector2<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;

    return *this;
  };

  // equality operator
  // -----------------

  bool operator== ( const vector2<X> &B )
  {
    for ( size_t i = 0;  i < 2 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class vector2

// arithmetic operators: vector2 = vector2 ? vector2
  // --------------------------------------------

template <class X> vector2<X> operator* (const vector2<X> &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

template <class X> vector2<X> operator/ (const vector2<X> &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];

  return C;
}

template <class X> vector2<X> operator+ (const vector2<X> &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

template <class X> vector2<X> operator- (const vector2<X> &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// arithmetic operators: vector2 = vector2 ? scalar
  // --------------------------------------------

template <class X> vector2<X> operator* (const vector2<X> &A, const X &B)
{
  vector2<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

template <class X> vector2<X> operator/ (const vector2<X> &A, const X &B)
{
  vector2<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

template <class X> vector2<X> operator+ (const vector2<X> &A, const X &B)
{
  vector2<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;

  return C;
}

template <class X> vector2<X> operator- (const vector2<X> &A, const X &B)
{
  vector2<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;

  return C;
}

// arithmetic operators: vector2 = scalar ? vector2
  // --------------------------------------------

template <class X> vector2<X> operator* (const X &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

template <class X> vector2<X> operator/ (const X &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];

  return C;
}

template <class X> vector2<X> operator+ (const X &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];

  return C;
}

template <class X> vector2<X> operator- (const X &A, const vector2<X> &B)
{
  vector2<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];

  return C;
}

// =================================================================================================
// print to screen
// =================================================================================================

template<class X> void inline tensor2_4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(3).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, m_data[i*8+j*4+k*2+l] );

}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2_4<X>& src)
{
  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          out << "(" << i << "," << j << "," << k << "," << l << ") " << src(i,j,k,l) << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2_2<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1]);
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[2],m_data[3]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2_2<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2_2s<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1]);
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[1],m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2_2s<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2_2d<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[0],m_data[2]);
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[2],m_data[1]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2_2d<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline vector2<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector2<X>& src)
{
  out << src(0) << ", " << src(1) << ";" << std::endl;

  return out;
}

// =================================================================================================
// identity tensors
// =================================================================================================

tensor2_4<double> inline identity2_4()
{
  tensor2_4<double> I(0.0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          if ( i==l and j==k )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor2_4<double> inline identity2_4rt()
{
  tensor2_4<double> I(0.0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          if ( i==k and j==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor2_4<double> inline identity2_II()
{
  tensor2_4<double> I(0.0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          if ( i==j and k==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor2_4<double> inline identity2_4s()
{ return (identity2_4()+identity2_4rt())/2.; }

// -------------------------------------------------------------------------------------------------

tensor2_4<double> inline identity2_4d()
{ return identity2_4s()-identity2_II()/static_cast<double>(2); }

// -------------------------------------------------------------------------------------------------

tensor2_2d<double> inline identity2_2()
{
  tensor2_2d<double> I(1.);

  return I;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X> tensor2_4<X> inline tensor2_4<X>::ddot(const tensor2_4<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          for ( size_t m=0; m<2; ++m )
            for ( size_t n=0; n<2; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l)*B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_4<X>::ddot(const tensor2_2<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_4<X>::ddot(const tensor2_2s<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_4<X>::ddot(const tensor2_2d<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,j) += (*this)(i,j,k,k)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::ddot(const tensor2_4<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2<X>::ddot(const tensor2_2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2<X>::ddot(const tensor2_2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2<X>::ddot(const tensor2_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[3]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2s<X>::ddot(const tensor2_4<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2s<X>::ddot(const tensor2_2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2s<X>::ddot(const tensor2_2s<X> &B) const
{
  X C;

  C  = m_data[0] * B[0];
  C += m_data[1] * B[1] * static_cast<X>(2);
  C += m_data[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2s<X>::ddot(const tensor2_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[2]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2d<X>::ddot(const tensor2_4<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      for ( size_t l=0; l<2; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2d<X>::ddot(const tensor2_2<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[3];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2d<X>::ddot(const tensor2_2s<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2d<X>::ddot(const tensor2_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[1];
}


// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::dot(const tensor2_2<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::dot(const tensor2_2s<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::dot(const tensor2_2d<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline tensor2_2<X>::dot(const vector2<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2s<X>::dot(const tensor2_2<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2s<X>::dot(const tensor2_2s<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2s<X>::dot(const tensor2_2d<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline tensor2_2s<X>::dot(const vector2<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2d<X>::dot(const tensor2_2<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2d<X>::dot(const tensor2_2s<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2d<X> inline tensor2_2d<X>::dot(const tensor2_2d<X> &B) const
{
  tensor2_2d<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    C(i,i) += (*this)(i,i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline tensor2_2d<X>::dot(const vector2<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    C(i) += (*this)(i,i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline vector2<X>::dot(const tensor2_2<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline vector2<X>::dot(const tensor2_2s<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector2<X> inline vector2<X>::dot(const tensor2_2d<X> &B) const
{
  vector2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    C(i) += (*this)(i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline vector2<X>::dot(const vector2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<2; ++i )
    C += (*this)(i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2<X>::dyadic(const tensor2_2<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2<X>::dyadic(const tensor2_2s<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2<X>::dyadic(const tensor2_2d<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2s<X>::dyadic(const tensor2_2<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2s<X>::dyadic(const tensor2_2s<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2s<X>::dyadic(const tensor2_2d<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2d<X>::dyadic(const tensor2_2<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      for ( size_t l=0; l<2; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2d<X>::dyadic(const tensor2_2s<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      for ( size_t l=0; l<2; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_2d<X>::dyadic(const tensor2_2d<X> &B) const
{
  tensor2_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t k=0; k<2; ++k )
      C(i,i,k,k) += (*this)(i,i)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline vector2<X>::dyadic(const vector2<X> &B) const
{
  tensor2_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      C(i,j) += (*this)(i)*B(j);

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X> tensor2_4<X> inline tensor2_4<X>::T() const
{
  tensor2_4<X> C;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_4<X>::RT() const
{
  tensor2_4<X> C;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_4<X> inline tensor2_4<X>::LT() const
{
  tensor2_4<X> C;

  for ( size_t i=0; i<2; ++i )
    for ( size_t j=0; j<2; ++j )
      for ( size_t k=0; k<2; ++k )
        for ( size_t l=0; l<2; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::T() const
{
  tensor2_2<X> C;

  C[0] = m_data[0];
  C[2] = m_data[1];
  C[1] = m_data[2];
  C[3] = m_data[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2s<X> inline tensor2_2s<X>::T() const
{
  tensor2_2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2d<X> inline tensor2_2d<X>::T() const
{
  tensor2_2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];

  return C;
}

// =================================================================================================
// miscellaneous
// =================================================================================================

template<class X> X inline tensor2_2<X>::trace() const
{
  return m_data[0]+m_data[3];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2s<X>::trace() const
{
  return m_data[0]+m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2d<X>::trace() const
{
  return m_data[0]+m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2<X>::det() const
{
  return m_data[0] * m_data[3] - m_data[1] * m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2s<X>::det() const
{
  return m_data[0] * m_data[2] - m_data[1] * m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2_2d<X>::det() const
{
  return m_data[0] * m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2<X> inline tensor2_2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2_2<X> C;

  C[0] =                      m_data[3] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] = static_cast<X>(-1) * m_data[2] / D;
  C[3] =                      m_data[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2s<X> inline tensor2_2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2_2s<X> C;

  C[0] =                      m_data[2] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] =                      m_data[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2_2d<X> inline tensor2_2d<X>::inv() const
{
  // allocate result
  tensor2_2d<X> C;

  C[0] = static_cast<X>(1) / m_data[0];
  C[1] = static_cast<X>(1) / m_data[1];

  return C;
}

// =================================================================================================
// create aliases to call class functions as functions, not members
// =================================================================================================

// products

// --

template<class X> tensor2_4 <X> inline ddot  (const tensor2_4 <X> &A, const tensor2_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_4 <X> &A, const tensor2_2 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_4 <X> &A, const tensor2_2s<X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_4 <X> &A, const tensor2_2d<X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_2 <X> &A, const tensor2_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_2s<X> &A, const tensor2_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor2_2 <X> inline ddot  (const tensor2_2d<X> &A, const tensor2_4 <X> &B)
{ return A.ddot(B); }

// --

template<class X>            X  inline ddot  (const tensor2_2 <X> &A, const tensor2_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{ return A.ddot(B); }

// --

template<class X> tensor2_2 <X> inline dot   (const tensor2_2 <X> &A, const tensor2_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor2_2 <X> inline dot   (const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor2_2d<X> inline dot   (const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{ return A.dot(B); }

// --

template<class X> vector2   <X> inline dot   (const tensor2_2 <X> &A, const vector2   <X> &B)
{ return A.dot(B); }

template<class X> vector2   <X> inline dot   (const tensor2_2s<X> &A, const vector2   <X> &B)
{ return A.dot(B); }

template<class X> vector2   <X> inline dot   (const tensor2_2d<X> &A, const vector2   <X> &B)
{ return A.dot(B); }

template<class X> vector2   <X> inline dot   (const vector2   <X> &A, const tensor2_2 <X> &B)
{ return A.dot(B); }

template<class X> vector2   <X> inline dot   (const vector2   <X> &A, const tensor2_2s<X> &B)
{ return A.dot(B); }

template<class X> vector2   <X> inline dot   (const vector2   <X> &A, const tensor2_2d<X> &B)
{ return A.dot(B); }

template<class X>            X  inline dot   (const vector2   <X> &A, const vector2   <X> &B)
{ return A.dot(B); }

// --

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2 <X> &A, const tensor2_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2 <X> &A, const tensor2_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2 <X> &A, const tensor2_2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2s<X> &A, const tensor2_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2s<X> &A, const tensor2_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2s<X> &A, const tensor2_2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2d<X> &A, const tensor2_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2d<X> &A, const tensor2_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor2_4 <X> inline dyadic(const tensor2_2d<X> &A, const tensor2_2d<X> &B)
{ return A.dyadic(B); }

// --

template<class X> tensor2_2 <X> inline dyadic(const vector2  <X> &A, const vector2  <X> &B)
{ return A.dyadic(B); }

// operations
template<class X> tensor2_2 <X> inline transpose (const tensor2_2 <X> &A) { return A.T    (); }
template<class X> tensor2_2s<X> inline transpose (const tensor2_2s<X> &A) { return A.T    (); }
template<class X> tensor2_2d<X> inline transpose (const tensor2_2d<X> &A) { return A.T    (); }
template<class X> tensor2_4 <X> inline transpose (const tensor2_4 <X> &A) { return A.T    (); }
template<class X> tensor2_4 <X> inline transposeR(const tensor2_4 <X> &A) { return A.RT   (); }
template<class X> tensor2_4 <X> inline transposeL(const tensor2_4 <X> &A) { return A.LT   (); }
template<class X> tensor2_2 <X> inline inv       (const tensor2_2 <X> &A) { return A.inv  (); }
template<class X> tensor2_2s<X> inline inv       (const tensor2_2s<X> &A) { return A.inv  (); }
template<class X> tensor2_2d<X> inline inv       (const tensor2_2d<X> &A) { return A.inv  (); }
template<class X>             X  inline det      (const tensor2_2 <X> &A) { return A.det  (); }
template<class X>             X  inline det      (const tensor2_2s<X> &A) { return A.det  (); }
template<class X>             X  inline det      (const tensor2_2d<X> &A) { return A.det  (); }
template<class X>             X  inline trace    (const tensor2_2 <X> &A) { return A.trace(); }
template<class X>             X  inline trace    (const tensor2_2s<X> &A) { return A.trace(); }
template<class X>             X  inline trace    (const tensor2_2d<X> &A) { return A.trace(); }

// =================================================================================================

} // namespace tensor

#endif

