/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef TENSOR3_H
#define TENSOR3_H

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

template<class X> class tensor3_4;
template<class X> class tensor3_2;
template<class X> class tensor3_2s;
template<class X> class tensor3_2d;
template<class X> class vector3;

// =================================================================================================
// cppmat::tensor3_4
// =================================================================================================

template<class X> class tensor3_4
{
private:

  X m_data[81]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor3_4(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor3_4(      X  D) { for ( size_t i=0; i<81; ++i ) m_data[i]=D;    };
  tensor3_4(const X *D) { for ( size_t i=0; i<81; ++i ) m_data[i]=D[i]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(4,3,bytes); };

  // copy constructor
  // ----------------

  // copy "tensor3_4" -> "tensor3_4" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_4<U> () const
  {
    tensor3_4<U> out;

    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
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

  size_t size() const { return 81; };
  size_t ndim() const { return 3;  };

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      C += std::abs(m_data[i  ]);
      C += std::abs(m_data[i+1]);
      C += std::abs(m_data[i+2]);
    }

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(0); }
  void ones        (     ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(1); }
  void setConstant ( X D ) { for ( size_t i=0; i<81; ++i ) m_data[i] = D;                 }

  // tensor products / operations
  // ----------------------------

  tensor3_4 <X> inline ddot(const tensor3_4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor3_2 <X> inline ddot(const tensor3_2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor3_2 <X> inline ddot(const tensor3_2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor3_2 <X> inline ddot(const tensor3_2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor3_4 <X> inline T   (                      ) const; // transposition   : B_lkji = A_ijkl
  tensor3_4 <X> inline RT  (                      ) const; // transposition   : B_ijlk = A_ijkl
  tensor3_4 <X> inline LT  (                      ) const; // transposition   : B_jikl = A_ijkl

  // index operators
  // ---------------

  X& operator[](size_t i)
  { return m_data[i]; };

  const X& operator[](size_t i) const
  { return m_data[i]; };

  X& operator()(size_t i, size_t j, size_t k, size_t l)
  { return m_data[i*27+j*9+k*3+l]; };

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*27+j*9+k*3+l]; };

  // arithmetic operators: tensor3_4 ?= tensor3_4
  // --------------------------------------------

  tensor3_4<X>& operator*= (const tensor3_4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] *= B[i  ];
      m_data[i+1] *= B[i+1];
      m_data[i+2] *= B[i+2];
    }

    return *this;
  };

  tensor3_4<X>& operator/= (const tensor3_4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] /= B[i  ];
      m_data[i+1] /= B[i+1];
      m_data[i+2] /= B[i+2];
    }

    return *this;
  };

  tensor3_4<X>& operator+= (const tensor3_4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] += B[i  ];
      m_data[i+1] += B[i+1];
      m_data[i+2] += B[i+2];
    }

    return *this;
  };

  tensor3_4<X>& operator-= (const tensor3_4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] -= B[i  ];
      m_data[i+1] -= B[i+1];
      m_data[i+2] -= B[i+2];
    }

    return *this;
  };

  // arithmetic operators: tensor3_4 ?= tensor3_4
  // --------------------------------------------

  tensor3_4<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] *= B;
      m_data[i+1] *= B;
      m_data[i+2] *= B;
    }

    return *this;
  };

  tensor3_4<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] /= B;
      m_data[i+1] /= B;
      m_data[i+2] /= B;
    }

    return *this;
  };

  tensor3_4<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] += B;
      m_data[i+1] += B;
      m_data[i+2] += B;
    }

    return *this;
  };

  tensor3_4<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] -= B;
      m_data[i+1] -= B;
      m_data[i+2] -= B;
    }

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor3_4<X> &B )
  {
    for ( size_t i = 0 ; i<81 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class tensor3_4

// arithmetic operators: tensor3_4 = tensor3_4 ? tensor3_4
// -------------------------------------------------------

template <class X> tensor3_4<X> operator* (const tensor3_4<X> &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator/ (const tensor3_4<X> &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator+ (const tensor3_4<X> &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator- (const tensor3_4<X> &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// arithmetic operators: tensor3_4 = tensor3_4 ? scalar
// ----------------------------------------------------

template <class X> tensor3_4<X> operator* (const tensor3_4<X> &A, const X &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C; }

template <class X> tensor3_4<X> operator/ (const tensor3_4<X> &A, const X &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C; }

template <class X> tensor3_4<X> operator+ (const tensor3_4<X> &A, const X &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C; }

template <class X> tensor3_4<X> operator- (const tensor3_4<X> &A, const X &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C; }

// arithmetic operators: tensor3_4 = scalar ? tensor3_4
// ----------------------------------------------------

template <class X> tensor3_4<X> operator* (const X &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator/ (const X &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator+ (const X &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

template <class X> tensor3_4<X> operator- (const X &A, const tensor3_4<X> &B)
{
  tensor3_4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// =================================================================================================
// cppmat::tensor3_2
// =================================================================================================

template<class X> class tensor3_2
{
private:

  X m_data[9]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor3_2(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  tensor3_2(      X  D) { for ( size_t i=0; i<9; ++i ) m_data[i]=D;    };
  tensor3_2(const X *D) { for ( size_t i=0; i<9; ++i ) m_data[i]=D[i]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,3,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor3_2" -> "tensor3_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2<U> () const
  {
    tensor3_2<U> out;

    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
    }

    return out;
  }

  // convert "tensor3_2 -> tensor3_2s"
  // WARNING: the output is symmetrized: "out(i,j) = ( this(i,j) + this(j,i) ) / 2."
  tensor3_2s<X> astensor2s()
  {
    tensor3_2s<X> out;

    out[0] =   m_data[0];
    out[1] = ( m_data[1] + m_data[3] ) / static_cast<X>(2);
    out[2] = ( m_data[2] + m_data[6] ) / static_cast<X>(2);
    out[3] =   m_data[4];
    out[4] = ( m_data[5] + m_data[7] ) / static_cast<X>(2);
    out[5] =   m_data[8];

    return out;
  }

  // convert "tensor3_2 -> tensor3_2d"
  // WARNING: all off-diagonal are discarded
  tensor3_2d<X> astensor2d()
  {
    tensor3_2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[4];
    out[2] = m_data[8];

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

  size_t size() const { return 9; };
  size_t ndim() const { return 3; };

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      C += std::abs(m_data[i  ]);
      C += std::abs(m_data[i+1]);
      C += std::abs(m_data[i+2]);
    }

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros       (     ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(0); };
  void ones        (     ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(1); };
  void setConstant ( X D ) { for ( size_t i=0; i<9; ++i ) m_data[i] = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor3_2 <X> inline dot   (const tensor3_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2 <X> inline dot   (const tensor3_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2 <X> inline dot   (const tensor3_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector3   <X> inline dot   (const vector3   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor3_2 <X> inline ddot  (const tensor3_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor3_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor3_4 <X> inline dyadic(const tensor3_2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_2 <X> inline T     (                      ) const; // transpose       : B_ij   = A_ji
  X             inline trace (                      ) const; // trace           : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor3_2 <X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i          )       { return m_data[i];     };
  const X& operator[](size_t i          ) const { return m_data[i];     };
  X&       operator()(size_t i, size_t j)       { return m_data[i*3+j]; };
  const X& operator()(size_t i, size_t j) const { return m_data[i*3+j]; };

  // arithmetic operators: tensor3_2 ?= tensor3_2
  // --------------------------------------------

  tensor3_2<X>& operator*= (const tensor3_2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] *= B[i  ];
      m_data[i+1] *= B[i+1];
      m_data[i+2] *= B[i+2];
    }

    return *this;
  };

  tensor3_2<X>& operator/= (const tensor3_2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] /= B[i  ];
      m_data[i+1] /= B[i+1];
      m_data[i+2] /= B[i+2];
    }

    return *this;
  };

  tensor3_2<X>& operator+= (const tensor3_2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] += B[i  ];
      m_data[i+1] += B[i+1];
      m_data[i+2] += B[i+2];
    }

    return *this;
  };

  tensor3_2<X>& operator-= (const tensor3_2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] -= B[i  ];
      m_data[i+1] -= B[i+1];
      m_data[i+2] -= B[i+2];
    }

    return *this;
  };

  // arithmetic operators: tensor3_2 ?= tensor3_2s
  // ---------------------------------------------

  tensor3_2<X>& operator*= (const tensor3_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1]; m_data[3] *= B[1];
    m_data[2] *= B[2]; m_data[6] *= B[2];
    m_data[4] *= B[3];
    m_data[5] *= B[4]; m_data[7] *= B[4];
    m_data[8] *= B[5];

    return *this;
  };

  tensor3_2<X>& operator/= (const tensor3_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1]; m_data[3] /= B[1];
    m_data[2] /= B[2]; m_data[6] /= B[2];
    m_data[4] /= B[3];
    m_data[5] /= B[4]; m_data[7] /= B[4];
    m_data[8] /= B[5];

    return *this;
  };

  tensor3_2<X>& operator+= (const tensor3_2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1]; m_data[3] += B[1];
    m_data[2] += B[2]; m_data[6] += B[2];
    m_data[4] += B[3];
    m_data[5] += B[4]; m_data[7] += B[4];
    m_data[8] += B[5];

    return *this;
  };

  tensor3_2<X>& operator-= (const tensor3_2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1]; m_data[3] -= B[1];
    m_data[2] -= B[2]; m_data[6] -= B[2];
    m_data[4] -= B[3];
    m_data[5] -= B[4]; m_data[7] -= B[4];
    m_data[8] -= B[5];

    return *this;
  };

  // arithmetic operators: tensor3_2 ?= tensor3_2d
  // -----------------------------------------

  tensor3_2<X>& operator*= (const tensor3_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[4] *= B[1];
    m_data[8] *= B[2];
    m_data[1] = m_data[2] = m_data[3] = m_data[5] = m_data[6] = m_data[7] = static_cast<X>(0);

    return *this;
  };

  tensor3_2<X>& operator+= (const tensor3_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[4] += B[1];
    m_data[8] += B[2];

    return *this;
  };

  tensor3_2<X>& operator-= (const tensor3_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[4] -= B[1];
    m_data[8] -= B[2];

    return *this;
  };

  // arithmetic operators: tensor3_2 ?= scalar
  // ---------------------------------------

  tensor3_2<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] *= B;
      m_data[i+1] *= B;
      m_data[i+2] *= B;
    }

    return *this;
  };

  tensor3_2<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] /= B;
      m_data[i+1] /= B;
      m_data[i+2] /= B;
    }

    return *this;
  };

  tensor3_2<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] += B;
      m_data[i+1] += B;
      m_data[i+2] += B;
    }

    return *this;
  };

  tensor3_2<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] -= B;
      m_data[i+1] -= B;
      m_data[i+2] -= B;
    }

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor3_2<X> &B )
  {
    for ( size_t i = 0;  i < 9 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor3_2s<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  };

  bool operator== ( const tensor3_2d<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool issymmetric()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i+1 ; j < 3 ; ++j )
        if ( m_data[ i*3 + j ] != m_data[ j*3 + i ] )
          return false;

    return true;
  };

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( i != j )
          if ( m_data[ i*3 + j ] )
            return false;

    return true;
  };

}; // class tensor3_2

// arithmetic operators: tensor3_2 = tensor3_2 ? tensor3_2
// -------------------------------------------------------

template <class X> tensor3_2<X> operator* (const tensor3_2<X> &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator/ (const tensor3_2<X> &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator+ (const tensor3_2<X> &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator- (const tensor3_2<X> &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// arithmetic operators: tensor3_2 = tensor3_2 ? tensor3_2s
// --------------------------------------------------------

template <class X> tensor3_2 <X> operator* (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[3] * B[1];
  C[2] = A[2] * B[2]; C[6] = A[6] * B[2];
  C[4] = A[4] * B[3];
  C[5] = A[5] * B[4]; C[7] = A[7] * B[4];
  C[8] = A[8] * B[5];

  return C;
}

template <class X> tensor3_2 <X> operator/ (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[3] / B[1];
  C[2] = A[2] / B[2]; C[6] = A[6] / B[2];
  C[4] = A[4] / B[3];
  C[5] = A[5] / B[4]; C[7] = A[7] / B[4];
  C[8] = A[8] / B[5];

  return C;
}

template <class X> tensor3_2 <X> operator+ (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[3] + B[1];
  C[2] = A[2] + B[2]; C[6] = A[6] + B[2];
  C[4] = A[4] + B[3];
  C[5] = A[5] + B[4]; C[7] = A[7] + B[4];
  C[8] = A[8] + B[5];

  return C;
}

template <class X> tensor3_2 <X> operator- (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[3] - B[1];
  C[2] = A[2] - B[2]; C[6] = A[6] - B[2];
  C[4] = A[4] - B[3];
  C[5] = A[5] - B[4]; C[7] = A[7] - B[4];
  C[8] = A[8] - B[5];

  return C;
}


// arithmetic operators: tensor3_2 = tensor3_2 ? tensor3_2d
// --------------------------------------------------------

template <class X> tensor3_2 <X> operator+ (const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3];
  C[4] = A[4] + B[1];
  C[5] = A[5];
  C[6] = A[6];
  C[7] = A[7];
  C[8] = A[8] + B[2];

  return C;
}

template <class X> tensor3_2 <X> operator- (const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3];
  C[4] = A[4] - B[1];
  C[5] = A[5];
  C[6] = A[6];
  C[7] = A[7];
  C[8] = A[8] - B[2];

  return C;
}

// arithmetic operators: tensor3_2 = tensor3_2 ? scalar
// ----------------------------------------------------

template <class X> tensor3_2<X> operator* (const tensor3_2<X> &A, const X &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C;
}

template <class X> tensor3_2<X> operator/ (const tensor3_2<X> &A, const X &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C;
}

template <class X> tensor3_2<X> operator+ (const tensor3_2<X> &A, const X &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C;
}

template <class X> tensor3_2<X> operator- (const tensor3_2<X> &A, const X &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C;
}

// arithmetic operators: tensor3_2 = tensor3_2s ? tensor3_2
// --------------------------------------------------

template <class X> tensor3_2 <X> operator* (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[1] * B[3];
  C[2] = A[2] * B[2]; C[6] = A[2] * B[6];
  C[4] = A[3] * B[4];
  C[5] = A[4] * B[5]; C[7] = A[4] * B[7];
  C[8] = A[5] * B[8];

  return C;
}

template <class X> tensor3_2 <X> operator/ (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[1] / B[3];
  C[2] = A[2] / B[2]; C[6] = A[2] / B[6];
  C[4] = A[3] / B[4];
  C[5] = A[4] / B[5]; C[7] = A[4] / B[7];
  C[8] = A[5] / B[8];

  return C;
}

template <class X> tensor3_2 <X> operator+ (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[1] + B[3];
  C[2] = A[2] + B[2]; C[6] = A[2] + B[6];
  C[4] = A[3] + B[4];
  C[5] = A[4] + B[5]; C[7] = A[4] + B[7];
  C[8] = A[5] + B[8];

  return C;
}

template <class X> tensor3_2 <X> operator- (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[1] - B[3];
  C[2] = A[2] - B[2]; C[6] = A[2] - B[6];
  C[4] = A[3] - B[4];
  C[5] = A[4] - B[5]; C[7] = A[4] - B[7];
  C[8] = A[5] - B[8];

  return C;
}

// arithmetic operators: tensor3_2 = tensor3_2d ? tensor3_2
// --------------------------------------------------

template <class X> tensor3_2 <X> operator+ (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] =        B[3];
  C[4] = A[1] + B[4];
  C[5] =        B[5];
  C[6] =        B[6];
  C[7] =        B[7];
  C[8] = A[2] + B[8];

  return C;
}

template <class X> tensor3_2 <X> operator- (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2 <X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] =      - B[3];
  C[4] = A[1] - B[4];
  C[5] =      - B[5];
  C[6] =      - B[6];
  C[7] =      - B[7];
  C[8] = A[2] - B[8];

  return C;
}

// arithmetic operators: tensor3_2 = scalar ? tensor3_2
// ------------------------------------------------

template <class X> tensor3_2<X> operator* (const X &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator/ (const X &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator+ (const X &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

template <class X> tensor3_2<X> operator- (const X &A, const tensor3_2<X> &B)
{
  tensor3_2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// =================================================================================================
// cppmat::tensor3_2s (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor3_2s
{
private:

  X m_data[6]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor3_2s(){};

  // explicit constructor: set to constant "D"
  tensor3_2s(X D) { for ( size_t i=0; i<6; ++i ) m_data[i]=D; };

  // explicit constructor: from full matrix
  // WARNING: the input is symmetrized: "this(i,j) = ( in(i,j) + in(j,i) ) / 2."
  tensor3_2s(const X *D)
  {
    // check for symmetry (code eliminated if "NDEBUG" is defined at the beginning of the code)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = i+1 ; j < 3 ; ++j )
          assert( D[ i*3 + j ] == D[ j*3 + i ] );
    #endif

    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        m_data[ i*3 - (i-1)*i/2 + j - i ] = D[ i*3 + j ];
  };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor3_2s", but of "cppmat::tensor3_2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,3,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor3_2s" -> "tensor3_2s" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2s<U> () const
  {
    tensor3_2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );
    out[3] = static_cast<U>( m_data[3] );
    out[4] = static_cast<U>( m_data[4] );
    out[5] = static_cast<U>( m_data[5] );

    return out;
  }

  // copy "const tensor3_2s" -> "tensor3_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2<U> () const
  {
    tensor3_2<U> out;

    out[0]          = m_data[0];
    out[1] = out[3] = m_data[1];
    out[2] = out[6] = m_data[2];
    out[4]          = m_data[3];
    out[5] = out[7] = m_data[4];
    out[8]          = m_data[5];

    return out;
  }

  // convert "tensor3_2s -> tensor3_2d"
  // WARNING: all off-diagonal are discarded
  tensor3_2d<X> astensor2d()
  {
    tensor3_2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[3];
    out[2] = m_data[5];

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

  size_t size() const { return 6; };
  size_t ndim() const { return 3; };

  // norm
  // ----

  X norm() const
  {
    X C;

    C  = std::abs(m_data[0]);
    C += std::abs(m_data[1]);
    C += std::abs(m_data[2]);
    C += std::abs(m_data[3]);
    C += std::abs(m_data[4]);
    C += std::abs(m_data[5]);

    return C;
  }

  // initialize to zero/one
  // ----------------------

  void zeros() { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(0); };
  void ones () { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=D; };

  // tensor products / operations
  // ----------------------------

  tensor3_2 <X> inline dot   (const tensor3_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2 <X> inline dot   (const tensor3_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2 <X> inline dot   (const tensor3_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector3   <X> inline dot   (const vector3   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor3_2 <X> inline ddot  (const tensor3_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor3_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor3_4 <X> inline dyadic(const tensor3_2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor3_2s<X> inline T     (                      ) const; // transpose       : B_ij   = A_ji
  X             inline trace (                      ) const; // trace           : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor3_2s<X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if ( i == 0 ) {
      if ( j == 0 ) return m_data[0];
      if ( j == 1 ) return m_data[1];
      else          return m_data[2];
    }
    if ( i == 1 ) {
      if ( j == 0 ) return m_data[1];
      if ( j == 1 ) return m_data[3];
      else          return m_data[4];
    }
    else {
      if ( j == 0 ) return m_data[2];
      if ( j == 1 ) return m_data[4];
      else          return m_data[5];
    }
  }

  const X& operator()(size_t i, size_t j) const
  {
    if ( i == 0 ) {
      if ( j == 0 ) return m_data[0];
      if ( j == 1 ) return m_data[1];
      else          return m_data[2];
    }
    if ( i == 1 ) {
      if ( j == 0 ) return m_data[1];
      if ( j == 1 ) return m_data[3];
      else          return m_data[4];
    }
    else {
      if ( j == 0 ) return m_data[2];
      if ( j == 1 ) return m_data[4];
      else          return m_data[5];
    }
  }

  // arithmetic operators: tensor3_2s ?= tensor3_2s
  // ------------------------------------------

  tensor3_2s<X>& operator*= (const tensor3_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];
    m_data[3] *= B[3];
    m_data[4] *= B[4];
    m_data[5] *= B[5];

    return *this;
  };

  tensor3_2s<X>& operator/= (const tensor3_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];
    m_data[3] /= B[3];
    m_data[4] /= B[4];
    m_data[5] /= B[5];

    return *this;
  };

  tensor3_2s<X>& operator+= (const tensor3_2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];
    m_data[3] += B[3];
    m_data[4] += B[4];
    m_data[5] += B[5];

    return *this;
  };

  tensor3_2s<X>& operator-= (const tensor3_2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];
    m_data[3] -= B[3];
    m_data[4] -= B[4];
    m_data[5] -= B[5];

    return *this;
  };

  // arithmetic operators: tensor3_2s ?= tensor3_2d
  // ------------------------------------------

  tensor3_2s<X>& operator*= (const tensor3_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[3] *= B[1];
    m_data[5] *= B[2];
    m_data[1] = m_data[2] = m_data[4] = static_cast<X>(0);

    return *this;
  };

  tensor3_2s<X>& operator+= (const tensor3_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[3] += B[1];
    m_data[5] += B[2];

    return *this;
  };

  tensor3_2s<X>& operator-= (const tensor3_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[3] -= B[1];
    m_data[5] -= B[2];

    return *this;
  };

  // arithmetic operators: tensor3_2s ?= scalar
  // ----------------------------------------

  tensor3_2s<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;
    m_data[3] *= B;
    m_data[4] *= B;
    m_data[5] *= B;

    return *this;
  };

  tensor3_2s<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;
    m_data[3] /= B;
    m_data[4] /= B;
    m_data[5] /= B;

    return *this;
  };

  tensor3_2s<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;
    m_data[3] += B;
    m_data[4] += B;
    m_data[5] += B;

    return *this;
  };

  tensor3_2s<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;
    m_data[3] -= B;
    m_data[4] -= B;
    m_data[5] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor3_2s<X> &B )
  {
    for ( size_t i = 0; i < 6 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor3_2<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = i ; j < 3 ; ++j ) {
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(j,i) ) return false;
      }
    }

    return true;
  };

  bool operator== ( const tensor3_2d<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;

    return true;
  };

  // structure-check operators
  // -------------------------

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i+1 ; j < 3 ; ++j )
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] )
          return false;

    return true;
  };

}; // class tensor3_2s

// arithmetic operators: tensor3_2s = tensor3_2s ? tensor3_2s
// ----------------------------------------------------------

template <class X> tensor3_2s<X> operator* (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];
  C[4] = A[4] * B[4];
  C[5] = A[5] * B[5];

  return C;
}

template <class X> tensor3_2s<X> operator/ (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];
  C[4] = A[4] / B[4];
  C[5] = A[5] / B[5];

  return C;
}

template <class X> tensor3_2s<X> operator+ (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];
  C[4] = A[4] + B[4];
  C[5] = A[5] + B[5];

  return C;
}

template <class X> tensor3_2s<X> operator- (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];
  C[4] = A[4] - B[4];
  C[5] = A[5] - B[5];

  return C;
}

// arithmetic operators: tensor3_2s = tensor3_2s ? tensor3_2d
// ----------------------------------------------------------

template <class X> tensor3_2s<X> operator+ (const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0]+B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3]+B[1];
  C[4] = A[4];
  C[5] = A[5]+B[2];

  return C;
}

template <class X> tensor3_2s<X> operator- (const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0]-B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3]-B[1];
  C[4] = A[4];
  C[5] = A[5]-B[2];

  return C;
}

// arithmetic operators: tensor3_2s = tensor3_2s ? scalar
// ------------------------------------------------------

template <class X> tensor3_2s<X> operator* (const tensor3_2s<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;
  C[4] = A[4] * B;
  C[5] = A[5] * B;

  return C;
}

template <class X> tensor3_2s<X> operator/ (const tensor3_2s<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;
  C[4] = A[4] / B;
  C[5] = A[5] / B;

  return C;
}

template <class X> tensor3_2s<X> operator+ (const tensor3_2s<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;
  C[4] = A[4] + B;
  C[5] = A[5] + B;

  return C;
}

template <class X> tensor3_2s<X> operator- (const tensor3_2s<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;
  C[4] = A[4] - B;
  C[5] = A[5] - B;

  return C;
}

// arithmetic operators: tensor3_2s = tensor3_2d ? scalar
// ------------------------------------------------------

template <class X> tensor3_2s<X> operator+ (const tensor3_2d<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] =        B;
  C[3] = A[1] + B;
  C[4] =        B;
  C[5] = A[2] + B;

  return C;
}

template <class X> tensor3_2s<X> operator- (const tensor3_2d<X> &A, const X &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] =      - B;
  C[3] = A[1] - B;
  C[4] =      - B;
  C[5] = A[2] - B;

  return C;
}

// arithmetic operators: tensor3_2s = tensor3_2d ? tensor3_2s
// ----------------------------------------------------------

template <class X> tensor3_2s<X> operator+ (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];
  C[4] =        B[4];
  C[5] = A[2] + B[5];

  return C;
}

template <class X> tensor3_2s<X> operator- (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];
  C[4] =      - B[4];
  C[5] = A[2] - B[5];

  return C;
}

// arithmetic operators: tensor3_2s = scalar ? tensor3_2s
// ------------------------------------------------------

template <class X> tensor3_2s<X> operator* (const X &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];
  C[4] = A * B[4];
  C[5] = A * B[5];

  return C;
}

template <class X> tensor3_2s<X> operator/ (const X &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];
  C[4] = A / B[4];
  C[5] = A / B[5];

  return C;
}

template <class X> tensor3_2s<X> operator+ (const X &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];
  C[4] = A + B[4];
  C[5] = A + B[5];

  return C;
}

template <class X> tensor3_2s<X> operator- (const X &A, const tensor3_2s<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];
  C[4] = A - B[4];
  C[5] = A - B[5];

  return C;
}

// arithmetic operators: tensor3_2s = scalar ? tensor3_2d
// ------------------------------------------------------

template <class X> tensor3_2s<X> operator+ (const X &A, const tensor3_2d<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A + B[1];
  C[4] = A;
  C[5] = A + B[2];

  return C;
}

template <class X> tensor3_2s<X> operator- (const X &A, const tensor3_2d<X> &B)
{
  tensor3_2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A - B[1];
  C[4] = A;
  C[5] = A - B[2];

  return C;
}

// =================================================================================================
// cppmat::tensor3_2d (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor3_2d
{
private:

  X m_data[4]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  tensor3_2d(){ m_data[3] = static_cast<X>(0); };

  // explicit constructor: set to constant "D"
  tensor3_2d(X D) { m_data[0] = m_data[1] = m_data[2] = D; };

  // explicit constructor: from full matrix
  // WARNING: all off-diagonal are discarded
  tensor3_2d(const X *D)
  {
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = 0 ; j < 3 ; ++j )
          if ( i != j )
            assert( ! D[ i*3 + j ] );
    #endif

    for ( size_t i=0; i<3; ++i )
      m_data[ i ] = D[ i*3 + i ];
  };

  // return strides array (see above)
  // WARNING: strides do not coincide with storage of "cppmat::tensor3_2d", but of "cppmat::tensor3_2"
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(2,3,bytes); };

  // copy constructors
  // -----------------

  // copy "tensor3_2d" -> "tensor3_2d" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2d<U> () const
  {
    tensor3_2d<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // copy "const tensor3_2d" -> "tensor3_2" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2<U> () const
  {
    tensor3_2<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[4] = static_cast<U>( m_data[1] );
    out[8] = static_cast<U>( m_data[2] );

    out[1] = out[2] = out[3] = out[5] = out[6] = out[7] = static_cast<U>(0);

    return out;
  }

  // copy "const tensor3_2d" -> "tensor3_2s" ( + change of type )
  // WARNING: all off-diagonal are discarded
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator tensor3_2s<U> () const
  {
    tensor3_2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[3] = static_cast<U>( m_data[1] );
    out[5] = static_cast<U>( m_data[2] );

    out[1] = out[2] = out[4] = static_cast<U>(0);

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
  size_t ndim() const { return 3; };

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

  void zeros       (     ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); };
  void ones        (     ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0] = m_data[1] = m_data[2] = D;                 };

  // tensor products / operations
  // ----------------------------

  tensor3_2 <X> inline dot   (const tensor3_2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2 <X> inline dot   (const tensor3_2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor3_2d<X> inline dot   (const tensor3_2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector3   <X> inline dot   (const vector3   <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor3_2 <X> inline ddot  (const tensor3_4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X             inline ddot  (const tensor3_2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X             inline ddot  (const tensor3_2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor3_4 <X> inline dyadic(const tensor3_2 <X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2s<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor3_4 <X> inline dyadic(const tensor3_2d<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  tensor3_2d<X> inline T     (                      ) const; // transpose      : B_ij   = A_ji
  X             inline trace (                      ) const; // trace          : A_ii
  X             inline det   (                      ) const; // determinant (only in 2D/3D)
  tensor3_2d<X> inline inv   (                      ) const; // inverse     (only in 2D/3D)

  // index operators
  // ---------------

  X&       operator[](size_t i )       { return m_data[i]; };
  const X& operator[](size_t i ) const { return m_data[i]; };

  X&       operator()(size_t i, size_t j)
  {
    if (i == j) return m_data[i];
    else        return m_data[3];
  }

  const X& operator()(size_t i, size_t j) const
  {
    if (i == j) return m_data[i];
    else        return m_data[3];
  }

  // arithmetic operators: tensor3_2d ?= tensor3_2d
  // ----------------------------------------------

  tensor3_2d<X>& operator*= (const tensor3_2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];

    return *this;
  };

  tensor3_2d<X>& operator+= (const tensor3_2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];

    return *this;
  };

  tensor3_2d<X>& operator-= (const tensor3_2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];

    return *this;
  };

  // arithmetic operators: tensor3_2d ?= tensor3_2
  // ---------------------------------------------

  tensor3_2d<X>& operator*= (const tensor3_2 <X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[4];
    m_data[2] *= B[8];

    return *this;
  };

  tensor3_2d<X>& operator/= (const tensor3_2 <X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[4];
    m_data[2] /= B[8];

    return *this;
  };

  // arithmetic operators: tensor3_2d ?= tensor3_2s
  // ----------------------------------------------

  tensor3_2d<X>& operator*= (const tensor3_2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[3];
    m_data[2] *= B[5];

    return *this;
  };

  tensor3_2d<X>& operator/= (const tensor3_2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[3];
    m_data[2] /= B[5];

    return *this;
  };

  // arithmetic operators: tensor3_2d ?= scalar
  // ------------------------------------------

  tensor3_2d<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;

    return *this;
  };

  tensor3_2d<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;

    return *this;
  };

  tensor3_2d<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;

    return *this;
  };

  tensor3_2d<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;

    return *this;
  };

  // equality operators
  // ------------------

  bool operator== ( const tensor3_2d<X> &B )
  {
    for ( size_t i = 0;  i < 3 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

  bool operator== ( const tensor3_2<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = 0 ; j < 3 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

  bool operator== ( const tensor3_2s<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = i ; j < 3 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  };

}; // class tensor3_2d

// arithmetic operators: tensor3_2d = tensor3_2d ? tensor3_2d
// ----------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

template <class X> tensor3_2d<X> operator+ (const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

template <class X> tensor3_2d<X> operator- (const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// arithmetic operators: tensor3_2d = tensor3_2d ? tensor3_2
// ---------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[4];
  C[2] = A[2] * B[8];

  return C;
}

template <class X> tensor3_2d<X> operator/ (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[4];
  C[2] = A[2] / B[8];

  return C;
}

// arithmetic operators: tensor3_2d = tensor3_2d ? tensor3_2s
// ----------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];
  C[2] = A[2] * B[5];

  return C;
}

template <class X> tensor3_2d<X> operator/ (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];
  C[2] = A[2] / B[5];

  return C;
}


// arithmetic operators: tensor3_2d = tensor3_2d ? scalar
// ------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2d<X> &A, const X &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

template <class X> tensor3_2d<X> operator/ (const tensor3_2d<X> &A, const X &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// arithmetic operators: tensor3_2d = tensor3_2 ? tensor3_2d
// ---------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[4] * B[1];
  C[2] = A[8] * B[2];

  return C;
}

// arithmetic operators: tensor3_2d = tensor3_2s ? tensor3_2d
// ----------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];
  C[2] = A[5] * B[2];

  return C;
}


// arithmetic operators: tensor3_2d = scalar ? tensor3_2d
// ------------------------------------------------------

template <class X> tensor3_2d<X> operator* (const X &A, const tensor3_2d<X> &B)
{
  tensor3_2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// =================================================================================================
// cppmat::vector3
// =================================================================================================

template<class X> class vector3
{
private:

  X m_data[3]; // data array

public:

  // constructors
  // ------------

  // implicit constructor
  vector3(){};

  // explicit constructors: set to constant "D", or copy from array (specified as pointer "*D")
  vector3(      X  D) { m_data[0] = m_data[1] = m_data[2] = D; };
  vector3(const X *D) { m_data[0] = D[0]; m_data[1] = D[1]; m_data[2] = D[2]; };

  // return strides array (see above)
  std::vector<size_t> strides(bool bytes=false) const { return _strides<X>(1,3,bytes); };

  // copy constructor
  // ----------------

  // copy "vector3" -> "vector3" ( + change of type )
  template<typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type>
  operator vector3<U> () const
  {
    vector3<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

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
  size_t ndim() const { return 3; };

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

  void zeros       (     ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); };
  void ones        (     ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); };
  void setConstant ( X D ) { m_data[0] = m_data[1] = m_data[2] = D;                 };

  // tensor products / operations
  // ----------------------------

  X            inline dot   (const vector3   <X> &B) const; // dot    product: C   = A_i*B_i
  vector3  <X> inline dot   (const tensor3_2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector3  <X> inline dot   (const tensor3_2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector3  <X> inline dot   (const tensor3_2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor3_2<X> inline dyadic(const vector3   <X> &B) const; // dyadic product: C_ij = A_i*B_j
  vector3  <X> inline cross (const vector3   <X> &B) const; // cross product (only in 3D)

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; };
  const X& operator[](size_t i) const { return m_data[i]; };
  X&       operator()(size_t i)       { return m_data[i]; };
  const X& operator()(size_t i) const { return m_data[i]; };

  // arithmetic operators: vector3 ?= vector3
  // --------------------------------------

  vector3<X>& operator*= (const vector3<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];

    return *this;
  };

  vector3<X>& operator/= (const vector3<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];

    return *this;
  };

  vector3<X>& operator+= (const vector3<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];

    return *this;
  };

  vector3<X>& operator-= (const vector3<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];

    return *this;
  };

  // arithmetic operators: vector3 ?= scalar
  // --------------------------------------

  vector3<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;

    return *this;
  };

  vector3<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;

    return *this;
  };

  vector3<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;

    return *this;
  };

  vector3<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;

    return *this;
  };

  // equality operator
  // -----------------

  bool operator== ( const vector3<X> &B )
  {
    for ( size_t i = 0;  i < 3 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  };

}; // class vector3

// arithmetic operators: vector3 = vector3 ? vector3
  // --------------------------------------------

template <class X> vector3<X> operator* (const vector3<X> &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

template <class X> vector3<X> operator/ (const vector3<X> &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

template <class X> vector3<X> operator+ (const vector3<X> &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

template <class X> vector3<X> operator- (const vector3<X> &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// arithmetic operators: vector3 = vector3 ? scalar
  // --------------------------------------------

template <class X> vector3<X> operator* (const vector3<X> &A, const X &B)
{
  vector3<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

template <class X> vector3<X> operator/ (const vector3<X> &A, const X &B)
{
  vector3<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

template <class X> vector3<X> operator+ (const vector3<X> &A, const X &B)
{
  vector3<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

template <class X> vector3<X> operator- (const vector3<X> &A, const X &B)
{
  vector3<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// arithmetic operators: vector3 = scalar ? vector3
  // --------------------------------------------

template <class X> vector3<X> operator* (const X &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

template <class X> vector3<X> operator/ (const X &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

template <class X> vector3<X> operator+ (const X &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

template <class X> vector3<X> operator- (const X &A, const vector3<X> &B)
{
  vector3<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// =================================================================================================
// print to screen
// =================================================================================================

template<class X> void inline tensor3_4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(3).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, m_data[i*27+j*9+k*3+l] );

}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor3_4<X>& src)
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          out << "(" << i << "," << j << "," << k << "," << l << ") " << src(i,j,k,l) << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor3_2<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[3],m_data[4],m_data[5]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[6],m_data[7],m_data[8]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor3_2<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor3_2s<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[1],m_data[3],m_data[4]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[2],m_data[4],m_data[5]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor3_2s<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor3_2d<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[3],m_data[3]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[3],m_data[1],m_data[3]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[3],m_data[3],m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor3_2d<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline vector3<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector3<X>& src)
{
  out << src(0) << ", " << src(1) << ", " << src(2) << ";" << std::endl;

  return out;
}

// =================================================================================================
// identity tensors
// =================================================================================================

tensor3_4<double> inline identity3_4()
{
  tensor3_4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==l and j==k )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor3_4<double> inline identity3_4rt()
{
  tensor3_4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==k and j==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor3_4<double> inline identity3_II()
{
  tensor3_4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==j and k==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor3_4<double> inline identity3_4s()
{ return (identity3_4()+identity3_4rt())/2.; }

// -------------------------------------------------------------------------------------------------

tensor3_4<double> inline identity3_4d()
{ return identity3_4s()-identity3_II()/static_cast<double>(3); }

// -------------------------------------------------------------------------------------------------

tensor3_2d<double> inline identity3_2()
{
  tensor3_2d<double> I(1.);

  return I;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X> tensor3_4<X> inline tensor3_4<X>::ddot(const tensor3_4<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          for ( size_t m=0; m<3; ++m )
            for ( size_t n=0; n<3; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l)*B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_4<X>::ddot(const tensor3_2<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_4<X>::ddot(const tensor3_2s<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_4<X>::ddot(const tensor3_2d<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j) += (*this)(i,j,k,k)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::ddot(const tensor3_4<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2<X>::ddot(const tensor3_2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2<X>::ddot(const tensor3_2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2<X>::ddot(const tensor3_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[4]*B[1] + m_data[8]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2s<X>::ddot(const tensor3_4<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2s<X>::ddot(const tensor3_2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2s<X>::ddot(const tensor3_2s<X> &B) const
{
  X C;

  C  = m_data[0] * B[0];
  C += m_data[1] * B[1] * static_cast<X>(2);
  C += m_data[2] * B[2] * static_cast<X>(2);
  C += m_data[3] * B[3];
  C += m_data[4] * B[4] * static_cast<X>(2);
  C += m_data[5] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2s<X>::ddot(const tensor3_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[3]*B[1] + m_data[5]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2d<X>::ddot(const tensor3_4<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2d<X>::ddot(const tensor3_2<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[4] + m_data[2]*B[8];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2d<X>::ddot(const tensor3_2s<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[3] + m_data[2]*B[5];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2d<X>::ddot(const tensor3_2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[1] + m_data[2]*B[2];
}


// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::dot(const tensor3_2<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::dot(const tensor3_2s<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::dot(const tensor3_2d<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline tensor3_2<X>::dot(const vector3<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2s<X>::dot(const tensor3_2<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2s<X>::dot(const tensor3_2s<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2s<X>::dot(const tensor3_2d<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline tensor3_2s<X>::dot(const vector3<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2d<X>::dot(const tensor3_2<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2d<X>::dot(const tensor3_2s<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2d<X> inline tensor3_2d<X>::dot(const tensor3_2d<X> &B) const
{
  tensor3_2d<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i,i) += (*this)(i,i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline tensor3_2d<X>::dot(const vector3<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i) += (*this)(i,i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline vector3<X>::dot(const tensor3_2<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline vector3<X>::dot(const tensor3_2s<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline vector3<X>::dot(const tensor3_2d<X> &B) const
{
  vector3<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i) += (*this)(i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline vector3<X>::dot(const vector3<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    C += (*this)(i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2<X>::dyadic(const tensor3_2<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2<X>::dyadic(const tensor3_2s<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2<X>::dyadic(const tensor3_2d<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2s<X>::dyadic(const tensor3_2<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2s<X>::dyadic(const tensor3_2s<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2s<X>::dyadic(const tensor3_2d<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2d<X>::dyadic(const tensor3_2<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2d<X>::dyadic(const tensor3_2s<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_2d<X>::dyadic(const tensor3_2d<X> &B) const
{
  tensor3_4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,i,k,k) += (*this)(i,i)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline vector3<X>::dyadic(const vector3<X> &B) const
{
  tensor3_2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector3<X> inline vector3<X>::cross(const vector3<X> &B) const
{
  vector3<X> C;

  C[0] =                     m_data[1]*B[2]-B[1]*m_data[2] ;
  C[1] = static_cast<X>(-1)*(m_data[0]*B[2]-B[0]*m_data[2]);
  C[2] =                     m_data[0]*B[1]-B[0]*m_data[1] ;

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X> tensor3_4<X> inline tensor3_4<X>::T() const
{
  tensor3_4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_4<X>::RT() const
{
  tensor3_4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_4<X> inline tensor3_4<X>::LT() const
{
  tensor3_4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::T() const
{
  tensor3_2<X> C;

  C[0] = m_data[0];
  C[3] = m_data[1];
  C[6] = m_data[2];
  C[1] = m_data[3];
  C[4] = m_data[4];
  C[7] = m_data[5];
  C[2] = m_data[6];
  C[5] = m_data[7];
  C[8] = m_data[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2s<X> inline tensor3_2s<X>::T() const
{
  tensor3_2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];
  C[3] = m_data[3];
  C[4] = m_data[4];
  C[5] = m_data[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2d<X> inline tensor3_2d<X>::T() const
{
  tensor3_2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// =================================================================================================
// miscellaneous
// =================================================================================================

template<class X> X inline tensor3_2<X>::trace() const
{
  return m_data[0]+m_data[4]+m_data[8];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2s<X>::trace() const
{
  return m_data[0]+m_data[3]+m_data[5];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2d<X>::trace() const
{
  return m_data[0]+m_data[1]+m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2<X>::det() const
{
  return ( m_data[0] * m_data[4] * m_data[8] +
           m_data[1] * m_data[5] * m_data[6] +
           m_data[2] * m_data[3] * m_data[7] ) -
         ( m_data[2] * m_data[4] * m_data[6] +
           m_data[1] * m_data[3] * m_data[8] +
           m_data[0] * m_data[5] * m_data[7] );
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2s<X>::det() const
{
  return (                     m_data[0] * m_data[3] * m_data[5] +
           static_cast<X>(2) * m_data[1] * m_data[2] * m_data[4] ) -
         (                     m_data[4] * m_data[4] * m_data[0] +
                               m_data[2] * m_data[2] * m_data[3] +
                               m_data[1] * m_data[1] * m_data[5] );
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor3_2d<X>::det() const
{
  return m_data[0] * m_data[1] * m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2<X> inline tensor3_2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor3_2<X> C;

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

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2s<X> inline tensor3_2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor3_2s<X> C;

  C[0] = (m_data[3]*m_data[5]-m_data[4]*m_data[4]) / D;
  C[1] = (m_data[2]*m_data[4]-m_data[1]*m_data[5]) / D;
  C[2] = (m_data[1]*m_data[4]-m_data[2]*m_data[3]) / D;
  C[3] = (m_data[0]*m_data[5]-m_data[2]*m_data[2]) / D;
  C[4] = (m_data[2]*m_data[1]-m_data[0]*m_data[4]) / D;
  C[5] = (m_data[0]*m_data[3]-m_data[1]*m_data[1]) / D;
  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor3_2d<X> inline tensor3_2d<X>::inv() const
{
  // allocate result
  tensor3_2d<X> C;

  C[0] = static_cast<X>(1) / m_data[0];
  C[1] = static_cast<X>(1) / m_data[1];
  C[2] = static_cast<X>(1) / m_data[2];

  return C;
}

// =================================================================================================
// create aliases to call class functions as functions, not members
// =================================================================================================

// products

// --

template<class X> tensor3_4 <X> inline ddot  (const tensor3_4 <X> &A, const tensor3_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_4 <X> &A, const tensor3_2 <X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_4 <X> &A, const tensor3_2s<X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_4 <X> &A, const tensor3_2d<X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_2 <X> &A, const tensor3_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_2s<X> &A, const tensor3_4 <X> &B)
{ return A.ddot(B); }

template<class X> tensor3_2 <X> inline ddot  (const tensor3_2d<X> &A, const tensor3_4 <X> &B)
{ return A.ddot(B); }

// --

template<class X>            X  inline ddot  (const tensor3_2 <X> &A, const tensor3_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{ return A.ddot(B); }

template<class X>            X  inline ddot  (const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{ return A.ddot(B); }

// --

template<class X> tensor3_2 <X> inline dot   (const tensor3_2 <X> &A, const tensor3_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{ return A.dot(B); }

template<class X> tensor3_2 <X> inline dot   (const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{ return A.dot(B); }

template<class X> tensor3_2d<X> inline dot   (const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{ return A.dot(B); }

// --

template<class X> vector3   <X> inline dot   (const tensor3_2 <X> &A, const vector3   <X> &B)
{ return A.dot(B); }

template<class X> vector3   <X> inline dot   (const tensor3_2s<X> &A, const vector3   <X> &B)
{ return A.dot(B); }

template<class X> vector3   <X> inline dot   (const tensor3_2d<X> &A, const vector3   <X> &B)
{ return A.dot(B); }

template<class X> vector3   <X> inline dot   (const vector3   <X> &A, const tensor3_2 <X> &B)
{ return A.dot(B); }

template<class X> vector3   <X> inline dot   (const vector3   <X> &A, const tensor3_2s<X> &B)
{ return A.dot(B); }

template<class X> vector3   <X> inline dot   (const vector3   <X> &A, const tensor3_2d<X> &B)
{ return A.dot(B); }

template<class X>            X  inline dot   (const vector3   <X> &A, const vector3   <X> &B)
{ return A.dot(B); }

// --

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2 <X> &A, const tensor3_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2 <X> &A, const tensor3_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2 <X> &A, const tensor3_2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2s<X> &A, const tensor3_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2s<X> &A, const tensor3_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2s<X> &A, const tensor3_2d<X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2d<X> &A, const tensor3_2 <X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2d<X> &A, const tensor3_2s<X> &B)
{ return A.dyadic(B); }

template<class X> tensor3_4 <X> inline dyadic(const tensor3_2d<X> &A, const tensor3_2d<X> &B)
{ return A.dyadic(B); }

// --

template<class X> tensor3_2 <X> inline dyadic(const vector3  <X> &A, const vector3  <X> &B)
{ return A.dyadic(B); }

// --

template<class X> vector3   <X> inline cross (const vector3  <X> &A, const vector3  <X> &B)
{ return A.cross (B); }

// operations
template<class X> tensor3_2 <X> inline transpose (const tensor3_2 <X> &A) { return A.T    (); }
template<class X> tensor3_2s<X> inline transpose (const tensor3_2s<X> &A) { return A.T    (); }
template<class X> tensor3_2d<X> inline transpose (const tensor3_2d<X> &A) { return A.T    (); }
template<class X> tensor3_4 <X> inline transpose (const tensor3_4 <X> &A) { return A.T    (); }
template<class X> tensor3_4 <X> inline transposeR(const tensor3_4 <X> &A) { return A.RT   (); }
template<class X> tensor3_4 <X> inline transposeL(const tensor3_4 <X> &A) { return A.LT   (); }
template<class X> tensor3_2 <X> inline inv       (const tensor3_2 <X> &A) { return A.inv  (); }
template<class X> tensor3_2s<X> inline inv       (const tensor3_2s<X> &A) { return A.inv  (); }
template<class X> tensor3_2d<X> inline inv       (const tensor3_2d<X> &A) { return A.inv  (); }
template<class X>             X  inline det      (const tensor3_2 <X> &A) { return A.det  (); }
template<class X>             X  inline det      (const tensor3_2s<X> &A) { return A.det  (); }
template<class X>             X  inline det      (const tensor3_2d<X> &A) { return A.det  (); }
template<class X>             X  inline trace    (const tensor3_2 <X> &A) { return A.trace(); }
template<class X>             X  inline trace    (const tensor3_2s<X> &A) { return A.trace(); }
template<class X>             X  inline trace    (const tensor3_2d<X> &A) { return A.trace(); }

// =================================================================================================

} // namespace tensor

#endif

