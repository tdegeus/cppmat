/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR3_H
#define CPPMAT_TENSOR3_H

#include "tensor.h" // strides
#include "macros.h"

namespace cppmat {
namespace cartesian3d {

// =================================================================================================
// forward declaration
// =================================================================================================

template<class X> class tensor4;
template<class X> class tensor2;
template<class X> class tensor2s;
template<class X> class tensor2d;
template<class X> class vector;

// =================================================================================================
// cppmat::cartesian3d::tensor4
// =================================================================================================

template<class X> class tensor4
{
private:

  // local data storage
  // - local container with the data
  X m_container[81];
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data = &m_container[0];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is initialized
  tensor4(){}

  // allocate tensor, initialize to constant "D"
  tensor4(X D) { for ( size_t i=0; i<81; ++i ) m_data[i]=D; }

  // copy raw pointer, correct storage
  tensor4(const X *D)
  {
    for ( size_t i = 0 ; i < 81 ; ++i ) m_data[i] = D[i];
  }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(X *D) { m_data = D; }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(const X *D)
  {
    for ( size_t i = 0 ; i < 81 ; ++i ) m_data[i] = D[i];
  }

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor4" -> "tensor4" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor4<U> () const
  {
    tensor4<U> out;

    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
    }

    return out;
  }

  // copy constructors : copy to std::vector
  // ---------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> () const
  {
    std::vector<U> out(81);

    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
    }

    return out;
  }

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return 81; }
  size_t ndim() const { return 3;  }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(4,3);

    return out;
  }

  // storage strides, optionally in bytes
  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> out = { 27, 9, 3, 1 };

    out[0] *= sizeof(X);
    out[1] *= sizeof(X);
    out[2] *= sizeof(X);
    out[3] *= sizeof(X);

    return out;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];      }
  auto     begin() const { return &m_data[0];      }
  auto     end  () const { return &m_data[0] + 81; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

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

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<81; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<81; ++i ) m_data[i] = static_cast<X>(1); }

  // tensor products / operations
  // ----------------------------

  tensor4<X> inline ddot(const tensor4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2<X> inline ddot(const tensor2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> inline ddot(const tensor2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> inline ddot(const tensor2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor4<X> inline T   (                    ) const; // transposition   : B_lkji = A_ijkl
  tensor4<X> inline RT  (                    ) const; // transposition   : B_ijlk = A_ijkl
  tensor4<X> inline LT  (                    ) const; // transposition   : B_jikl = A_ijkl

  // index operators
  // ---------------

  X& operator[](size_t i)
  { return m_data[i]; }

  const X& operator[](size_t i) const
  { return m_data[i]; }

  X& operator()(size_t i, size_t j, size_t k, size_t l)
  { return m_data[i*27+j*9+k*3+l]; }

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*27+j*9+k*3+l]; }

  // arithmetic operators: tensor4 ?= tensor4
  // ----------------------------------------

  tensor4<X>& operator*= (const tensor4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] *= B[i  ];
      m_data[i+1] *= B[i+1];
      m_data[i+2] *= B[i+2];
    }

    return *this;
  }

  tensor4<X>& operator/= (const tensor4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] /= B[i  ];
      m_data[i+1] /= B[i+1];
      m_data[i+2] /= B[i+2];
    }

    return *this;
  }

  tensor4<X>& operator+= (const tensor4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] += B[i  ];
      m_data[i+1] += B[i+1];
      m_data[i+2] += B[i+2];
    }

    return *this;
  }

  tensor4<X>& operator-= (const tensor4<X> &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] -= B[i  ];
      m_data[i+1] -= B[i+1];
      m_data[i+2] -= B[i+2];
    }

    return *this;
  }

  // arithmetic operators: tensor4 ?= scalar
  // ---------------------------------------

  tensor4<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] *= B;
      m_data[i+1] *= B;
      m_data[i+2] *= B;
    }

    return *this;
  }

  tensor4<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] /= B;
      m_data[i+1] /= B;
      m_data[i+2] /= B;
    }

    return *this;
  }

  tensor4<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] += B;
      m_data[i+1] += B;
      m_data[i+2] += B;
    }

    return *this;
  }

  tensor4<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < 81 ; i += 3 ) {
      m_data[i  ] -= B;
      m_data[i+1] -= B;
      m_data[i+2] -= B;
    }

    return *this;
  }

  // equality operators: tensor4 == tensor4
  // --------------------------------------

  bool operator== ( const tensor4<X> &B )
  {
    for ( size_t i = 0 ; i<81 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

}; // class tensor4

// arithmetic operators: tensor4 = tensor4 ? tensor4
// -------------------------------------------------

template <class X> tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// arithmetic operators: tensor4 = tensor4 ? scalar
// ------------------------------------------------

template <class X> tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C; }

template <class X> tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C; }

template <class X> tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C; }

template <class X> tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C; }

// arithmetic operators: tensor4 = scalar ? tensor4
// ------------------------------------------------

template <class X> tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

template <class X> tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 ) {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// =================================================================================================
// cppmat::cartesian3d::tensor2
// =================================================================================================

template<class X> class tensor2
{
private:

  // local data storage
  // - local container with the data
  X m_container[9];
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data = &m_container[0];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is initialized
  tensor2(){}

  // allocate tensor, initialize to constant "D"
  tensor2(X D) { for ( size_t i=0; i<9; ++i ) m_data[i]=D; }

  // copy from raw pointer, correct storage
  tensor2(const X *D)
  {
    for ( size_t i = 0 ; i < 9 ; ++i )
      m_data[i] = D[i];
  }

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  tensor2(const Eigen::Matrix<X,3,3,Eigen::RowMajor> &D)
  {
    for ( size_t i = 0 ; i < 9 ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  tensor2(const Eigen::Matrix<X,3,3,Eigen::ColMajor> &D)
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        m_data[i*3+j] = D(i,j);
  }
  #endif

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(X *D) { m_data = D; }

  // pointer to Eigen array
  // N.B. only possible for matching 'RowMajor' storage, the user has to keep the pointer alive
  #ifdef CPPMAT_EIGEN
  void map(const Eigen::Matrix<X,3,3,Eigen::RowMajor> &D) { m_data = const_cast<X*>(&D(0)); }
  #endif

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(const X *D)
  {
    for ( size_t i = 0 ; i < 9 ; ++i )
      m_data[i] = D[i];
  }

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,3,3,Eigen::RowMajor> &D)
  {
    for ( size_t i = 0 ; i < 9 ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,3,3,Eigen::ColMajor> &D)
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        m_data[i*3+j] = D(i,j);
  }
  #endif

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out;

    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
    }

    return out;
  }

  // copy constructors : copy to std::vector
  // ---------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> () const
  {
    std::vector<U> out(9);

    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      out[i  ] = static_cast<U>( m_data[i  ] );
      out[i+1] = static_cast<U>( m_data[i+1] );
      out[i+2] = static_cast<U>( m_data[i+2] );
    }

    return out;
  }

  // copy constructors : copy to Eigen object
  // ----------------------------------------

  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,3,3,Eigen::RowMajor> () const
  {
    Eigen::Matrix<U,3,3,Eigen::RowMajor> out;

    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      out(i  ) = static_cast<U>( m_data[i  ] );
      out(i+1) = static_cast<U>( m_data[i+1] );
      out(i+2) = static_cast<U>( m_data[i+2] );
    }

    return out;
  }
  #endif

  // conversions : convert to incompatible tensor objects -> SOME TERMS ARE DISCARDED
  // --------------------------------------------------------------------------------

  // "tensor2 -> tensor2s" (output is symmetrized: "out(i,j) = ( this(i,j) + this(j,i) ) / 2.")
  tensor2s<X> astensor2s()
  {
    tensor2s<X> out;

    out[0] =   m_data[0];
    out[1] = ( m_data[1] + m_data[3] ) / static_cast<X>(2);
    out[2] = ( m_data[2] + m_data[6] ) / static_cast<X>(2);
    out[3] =   m_data[4];
    out[4] = ( m_data[5] + m_data[7] ) / static_cast<X>(2);
    out[5] =   m_data[8];

    return out;
  }

  // "tensor2 -> tensor2d" (all off-diagonal terms are discarded)
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[4];
    out[2] = m_data[8];

    return out;
  }

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return 9; }
  size_t ndim() const { return 4; }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(2,3);

    return out;
  }

  // storage strides, optionally in bytes
  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> out = { 3, 1 };

    out[0] *= sizeof(X);
    out[1] *= sizeof(X);

    return out;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];     }
  auto     begin() const { return &m_data[0];     }
  auto     end  () const { return &m_data[0] + 9; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

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

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<9; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<9; ++i ) m_data[i] = static_cast<X>(1); }

  // tensor products / operations
  // ----------------------------

  tensor2<X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2<X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2<X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2<X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X          inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X          inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X          inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4<X> inline dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4<X> inline dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4<X> inline dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2<X> inline T     (                    ) const; // transpose       : B_ij   = A_ji
  X          inline trace (                    ) const; // trace           : A_ii
  X          inline det   (                    ) const; // determinant
  tensor2<X> inline inv   (                    ) const; // inverse

  // index operators
  // ---------------

  X&       operator[](size_t i          )       { return m_data[i];     }
  const X& operator[](size_t i          ) const { return m_data[i];     }
  X&       operator()(size_t i, size_t j)       { return m_data[i*3+j]; }
  const X& operator()(size_t i, size_t j) const { return m_data[i*3+j]; }

  // arithmetic operators: tensor2 ?= tensor2
  // ----------------------------------------

  tensor2<X>& operator*= (const tensor2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] *= B[i  ];
      m_data[i+1] *= B[i+1];
      m_data[i+2] *= B[i+2];
    }

    return *this;
  }

  tensor2<X>& operator/= (const tensor2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] /= B[i  ];
      m_data[i+1] /= B[i+1];
      m_data[i+2] /= B[i+2];
    }

    return *this;
  }

  tensor2<X>& operator+= (const tensor2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] += B[i  ];
      m_data[i+1] += B[i+1];
      m_data[i+2] += B[i+2];
    }

    return *this;
  }

  tensor2<X>& operator-= (const tensor2<X> &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] -= B[i  ];
      m_data[i+1] -= B[i+1];
      m_data[i+2] -= B[i+2];
    }

    return *this;
  }

  // arithmetic operators: tensor2 ?= tensor2s
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1]; m_data[3] *= B[1];
    m_data[2] *= B[2]; m_data[6] *= B[2];
    m_data[4] *= B[3];
    m_data[5] *= B[4]; m_data[7] *= B[4];
    m_data[8] *= B[5];

    return *this;
  }

  tensor2<X>& operator/= (const tensor2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1]; m_data[3] /= B[1];
    m_data[2] /= B[2]; m_data[6] /= B[2];
    m_data[4] /= B[3];
    m_data[5] /= B[4]; m_data[7] /= B[4];
    m_data[8] /= B[5];

    return *this;
  }

  tensor2<X>& operator+= (const tensor2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1]; m_data[3] += B[1];
    m_data[2] += B[2]; m_data[6] += B[2];
    m_data[4] += B[3];
    m_data[5] += B[4]; m_data[7] += B[4];
    m_data[8] += B[5];

    return *this;
  }

  tensor2<X>& operator-= (const tensor2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1]; m_data[3] -= B[1];
    m_data[2] -= B[2]; m_data[6] -= B[2];
    m_data[4] -= B[3];
    m_data[5] -= B[4]; m_data[7] -= B[4];
    m_data[8] -= B[5];

    return *this;
  }

  // arithmetic operators: tensor2 ?= tensor2d
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[4] *= B[1];
    m_data[8] *= B[2];
    m_data[1] = m_data[2] = m_data[3] = m_data[5] = m_data[6] = m_data[7] = static_cast<X>(0);

    return *this;
  }

  tensor2<X>& operator+= (const tensor2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[4] += B[1];
    m_data[8] += B[2];

    return *this;
  }

  tensor2<X>& operator-= (const tensor2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[4] -= B[1];
    m_data[8] -= B[2];

    return *this;
  }

  // arithmetic operators: tensor2 ?= scalar
  // ---------------------------------------

  tensor2<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] *= B;
      m_data[i+1] *= B;
      m_data[i+2] *= B;
    }

    return *this;
  }

  tensor2<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] /= B;
      m_data[i+1] /= B;
      m_data[i+2] /= B;
    }

    return *this;
  }

  tensor2<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] += B;
      m_data[i+1] += B;
      m_data[i+2] += B;
    }

    return *this;
  }

  tensor2<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < 9 ; i += 3 ) {
      m_data[i  ] -= B;
      m_data[i+1] -= B;
      m_data[i+2] -= B;
    }

    return *this;
  }

  // equality operators: tensor2 == tensor2?
  // ---------------------------------------

  bool operator== ( const tensor2<X> &B )
  {
    for ( size_t i = 0;  i < 9 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2s<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  }

  bool operator== ( const tensor2d<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( m_data[i*3+j] != B(i,j) )
          return false;

    return true;
  }

  // structure-check operators
  // -------------------------

  bool issymmetric()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i+1 ; j < 3 ; ++j )
        if ( m_data[ i*3 + j ] != m_data[ j*3 + i ] )
          return false;

    return true;
  }

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = 0 ; j < 3 ; ++j )
        if ( i != j )
          if ( m_data[ i*3 + j ] )
            return false;

    return true;
  }

}; // class tensor2

// arithmetic operators: tensor2 = tensor2 ? tensor2
// -------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// arithmetic operators: tensor2 = tensor2 ? tensor2s
// --------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[3] * B[1];
  C[2] = A[2] * B[2]; C[6] = A[6] * B[2];
  C[4] = A[4] * B[3];
  C[5] = A[5] * B[4]; C[7] = A[7] * B[4];
  C[8] = A[8] * B[5];

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[3] / B[1];
  C[2] = A[2] / B[2]; C[6] = A[6] / B[2];
  C[4] = A[4] / B[3];
  C[5] = A[5] / B[4]; C[7] = A[7] / B[4];
  C[8] = A[8] / B[5];

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[3] + B[1];
  C[2] = A[2] + B[2]; C[6] = A[6] + B[2];
  C[4] = A[4] + B[3];
  C[5] = A[5] + B[4]; C[7] = A[7] + B[4];
  C[8] = A[8] + B[5];

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[3] - B[1];
  C[2] = A[2] - B[2]; C[6] = A[6] - B[2];
  C[4] = A[4] - B[3];
  C[5] = A[5] - B[4]; C[7] = A[7] - B[4];
  C[8] = A[8] - B[5];

  return C;
}


// arithmetic operators: tensor2 = tensor2 ? tensor2d
// --------------------------------------------------

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2<X> C;

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

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2 <X> C;

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

// arithmetic operators: tensor2 = tensor2 ? scalar
// ------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C;
}

// arithmetic operators: tensor2 = tensor2s ? tensor2
// --------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[1] * B[3];
  C[2] = A[2] * B[2]; C[6] = A[2] * B[6];
  C[4] = A[3] * B[4];
  C[5] = A[4] * B[5]; C[7] = A[4] * B[7];
  C[8] = A[5] * B[8];

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[1] / B[3];
  C[2] = A[2] / B[2]; C[6] = A[2] / B[6];
  C[4] = A[3] / B[4];
  C[5] = A[4] / B[5]; C[7] = A[4] / B[7];
  C[8] = A[5] / B[8];

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[1] + B[3];
  C[2] = A[2] + B[2]; C[6] = A[2] + B[6];
  C[4] = A[3] + B[4];
  C[5] = A[4] + B[5]; C[7] = A[4] + B[7];
  C[8] = A[5] + B[8];

  return C;
}

template <class X> tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[1] - B[3];
  C[2] = A[2] - B[2]; C[6] = A[2] - B[6];
  C[4] = A[3] - B[4];
  C[5] = A[4] - B[5]; C[7] = A[4] - B[7];
  C[8] = A[5] - B[8];

  return C;
}

// arithmetic operators: tensor2 = tensor2d ? tensor2
// --------------------------------------------------

template <class X> tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

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

template <class X> tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

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

// arithmetic operators: tensor2 = scalar ? tensor2
// ------------------------------------------------

template <class X> tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

template <class X> tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 ) {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// =================================================================================================
// cppmat::cartesian3d::tensor2s (symmetric storage of "cppmat::tensor3::tensor")
// =================================================================================================

template<class X> class tensor2s
{
private:

  // local data storage
  // - local container with the data
  X m_container[6];
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data = &m_container[0];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is initialized
  tensor2s(){}

  // allocate tensor, initialize to constant "D"
  tensor2s(X D) { for ( size_t i=0; i<6; ++i ) m_data[i]=D; }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(X *D) { m_data = D; }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(const X *D)
  {
    for ( size_t i = 0 ; i < 6 ; ++i )
      m_data[i] = D[i];
  }

  // raw pointer, stored as if it was "tensor2"
  void copyDense(const X *D)
  {
    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = i+1 ; j < 3 ; ++j )
          assert( D[i*3+j] == D[j*3+i] );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        m_data[i*3-(i-1)*i/3+j-i] = D[i*3+j];
  }

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,3,3,Eigen::RowMajor> &D)
  {
    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = i+1 ; j < 3 ; ++j )
          assert( D(i,j) == D(j,i) );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        m_data[i*3-(i-1)*i/3+j-i] = D(i,j);
  }
  #endif

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,3,3,Eigen::ColMajor> &D)
  {
    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = i+1 ; j < 3 ; ++j )
          assert( D(i,j) == D(j,i) );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        m_data[i*3-(i-1)*i/3+j-i] = D(i,j);
  }
  #endif

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2s" -> "tensor2s" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2s<U> () const
  {
    tensor2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );
    out[3] = static_cast<U>( m_data[3] );
    out[4] = static_cast<U>( m_data[4] );
    out[5] = static_cast<U>( m_data[5] );

    return out;
  }

  // copy "const tensor2s" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out;

    out[0]          = m_data[0];
    out[1] = out[3] = m_data[1];
    out[2] = out[6] = m_data[2];
    out[4]          = m_data[3];
    out[5] = out[7] = m_data[4];
    out[8]          = m_data[5];

    return out;
  }

  // copy constructors : copy to std::vector
  // ---------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> () const
  {
    std::vector<U> out(6);

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );
    out[3] = static_cast<U>( m_data[3] );
    out[4] = static_cast<U>( m_data[4] );
    out[5] = static_cast<U>( m_data[5] );

    return out;
  }

  // conversions : convert to incompatible tensor objects -> SOME TERMS ARE DISCARDED
  // --------------------------------------------------------------------------------

  // "tensor2s -> tensor2d" (all off-diagonal terms are discarded)
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out;

    out[0] = m_data[0];
    out[1] = m_data[3];
    out[2] = m_data[5];

    return out;
  }

  // dimensions
  // ----------

  size_t size() const { return 6; }
  size_t ndim() const { return 3; }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];     }
  auto     begin() const { return &m_data[0];     }
  auto     end  () const { return &m_data[0] + 6; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

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

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=D;                 }
  void setZero    (   ) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(0); }
  void setOnes    (   ) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(1); }
  void zeros      (   ) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(0); }
  void ones       (   ) { m_data[0]=m_data[1]=m_data[2]=m_data[3]=m_data[4]=m_data[5]=static_cast<X>(1); }

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
  X           inline det   (                    ) const; // determinant
  tensor2s<X> inline inv   (                    ) const; // inverse

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; }
  const X& operator[](size_t i) const { return m_data[i]; }

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

  // arithmetic operators: tensor2s ?= tensor2s
  // ------------------------------------------

  tensor2s<X>& operator*= (const tensor2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];
    m_data[3] *= B[3];
    m_data[4] *= B[4];
    m_data[5] *= B[5];

    return *this;
  }

  tensor2s<X>& operator/= (const tensor2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];
    m_data[3] /= B[3];
    m_data[4] /= B[4];
    m_data[5] /= B[5];

    return *this;
  }

  tensor2s<X>& operator+= (const tensor2s<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];
    m_data[3] += B[3];
    m_data[4] += B[4];
    m_data[5] += B[5];

    return *this;
  }

  tensor2s<X>& operator-= (const tensor2s<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];
    m_data[3] -= B[3];
    m_data[4] -= B[4];
    m_data[5] -= B[5];

    return *this;
  }

  // arithmetic operators: tensor2s ?= tensor2d
  // ------------------------------------------

  tensor2s<X>& operator*= (const tensor2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[3] *= B[1];
    m_data[5] *= B[2];
    m_data[1] = m_data[2] = m_data[4] = static_cast<X>(0);

    return *this;
  }

  tensor2s<X>& operator+= (const tensor2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[3] += B[1];
    m_data[5] += B[2];

    return *this;
  }

  tensor2s<X>& operator-= (const tensor2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[3] -= B[1];
    m_data[5] -= B[2];

    return *this;
  }

  // arithmetic operators: tensor2s ?= scalar
  // ----------------------------------------

  tensor2s<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;
    m_data[3] *= B;
    m_data[4] *= B;
    m_data[5] *= B;

    return *this;
  }

  tensor2s<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;
    m_data[3] /= B;
    m_data[4] /= B;
    m_data[5] /= B;

    return *this;
  }

  tensor2s<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;
    m_data[3] += B;
    m_data[4] += B;
    m_data[5] += B;

    return *this;
  }

  tensor2s<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;
    m_data[3] -= B;
    m_data[4] -= B;
    m_data[5] -= B;

    return *this;
  }

  // equality operators: tensor2s == tensor2?
  // ----------------------------------------

  bool operator== ( const tensor2s<X> &B )
  {
    for ( size_t i = 0; i < 6 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = i ; j < 3 ; ++j ) {
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(j,i) ) return false;
      }
    }

    return true;
  }

  bool operator== ( const tensor2d<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i ; j < 3 ; ++j )
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] != B(i,j) ) return false;

    return true;
  }

  // structure-check operators
  // -------------------------

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      for ( size_t j = i+1 ; j < 3 ; ++j )
        if ( m_data[ i*3 - (i-1)*i/2 + j - i ] )
          return false;

    return true;
  }

}; // class tensor2s

// arithmetic operators: tensor2s = tensor2s ? tensor2s
// ----------------------------------------------------

template <class X> tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];
  C[4] = A[4] * B[4];
  C[5] = A[5] * B[5];

  return C;
}

template <class X> tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];
  C[4] = A[4] / B[4];
  C[5] = A[5] / B[5];

  return C;
}

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];
  C[4] = A[4] + B[4];
  C[5] = A[5] + B[5];

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];
  C[4] = A[4] - B[4];
  C[5] = A[5] - B[5];

  return C;
}

// arithmetic operators: tensor2s = tensor2s ? tensor2d
// ----------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0]+B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3]+B[1];
  C[4] = A[4];
  C[5] = A[5]+B[2];

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0]-B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3]-B[1];
  C[4] = A[4];
  C[5] = A[5]-B[2];

  return C;
}

// arithmetic operators: tensor2s = tensor2s ? scalar
// --------------------------------------------------

template <class X> tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;
  C[4] = A[4] * B;
  C[5] = A[5] * B;

  return C;
}

template <class X> tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;
  C[4] = A[4] / B;
  C[5] = A[5] / B;

  return C;
}

template <class X> tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;
  C[4] = A[4] + B;
  C[5] = A[5] + B;

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;
  C[4] = A[4] - B;
  C[5] = A[5] - B;

  return C;
}

// arithmetic operators: tensor2s = tensor2d ? scalar
// --------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] =        B;
  C[3] = A[1] + B;
  C[4] =        B;
  C[5] = A[2] + B;

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] =      - B;
  C[3] = A[1] - B;
  C[4] =      - B;
  C[5] = A[2] - B;

  return C;
}

// arithmetic operators: tensor2s = tensor2d ? tensor2s
// ----------------------------------------------------

template <class X> tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];
  C[4] =        B[4];
  C[5] = A[2] + B[5];

  return C;
}

template <class X> tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];
  C[4] =      - B[4];
  C[5] = A[2] - B[5];

  return C;
}

// arithmetic operators: tensor2s = scalar ? tensor2s
// --------------------------------------------------

template <class X> tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];
  C[4] = A * B[4];
  C[5] = A * B[5];

  return C;
}

template <class X> tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];
  C[4] = A / B[4];
  C[5] = A / B[5];

  return C;
}

template <class X> tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];
  C[4] = A + B[4];
  C[5] = A + B[5];

  return C;
}

template <class X> tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];
  C[4] = A - B[4];
  C[5] = A - B[5];

  return C;
}

// arithmetic operators: tensor2s = scalar ? tensor2d
// --------------------------------------------------

template <class X> tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A + B[1];
  C[4] = A;
  C[5] = A + B[2];

  return C;
}

template <class X> tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A - B[1];
  C[4] = A;
  C[5] = A - B[2];

  return C;
}

// =================================================================================================
// cppmat::cartesian3d::tensor2d (symmetric storage of "cppmat::tensor3::tensor")
// =================================================================================================

template<class X> class tensor2d
{
private:

  // local data storage
  // - local container with the data (N.B. the final index is used as a dummy zero)
  X m_container[3];
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data = &m_container[0];

  // dummy parameter, used to return "0" for any off-diagonal entry
  X m_zero[1];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is initialized (and set the dummy parameter)
  tensor2d(){ m_zero[0] = static_cast<X>(0); }

  // allocate tensor, initialize to constant "D" (and set the dummy parameter)
  tensor2d(X D) { m_data[0] = m_data[1] = m_data[2] = D; m_zero[0] = static_cast<X>(0); }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(X *D) { m_data = D; }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(const X *D)
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      m_data[i] = D[i];
  }

  // raw pointer, stored as if it was "tensor2"
  void copyDense(const X *D)
  {
    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = 0 ; j < 3 ; ++j )
          if ( i != j )
            assert( ! D[i*3+j] );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < 3 ; ++i )
      m_data[i] = D[i*3+i];
  }

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,3,3,Eigen::RowMajor> &D)
  {
    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = 0 ; j < 3 ; ++j )
          if ( i != j )
            assert( ! D(i,j) );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < 3 ; ++i )
      m_data[i] = D(i,i);
  }
  #endif

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,3,3,Eigen::ColMajor> &D)
  {
    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < 3 ; ++i )
        for ( size_t j = 0 ; j < 3 ; ++j )
          if ( i != j )
            assert( ! D(i,j) );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < 3 ; ++i )
      m_data[i] = D(i,i);
  }
  #endif

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2d" -> "tensor2d" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2d<U> () const
  {
    tensor2d<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // copy "const tensor2d" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[4] = static_cast<U>( m_data[1] );
    out[8] = static_cast<U>( m_data[2] );

    out[1] = out[2] = out[3] = out[5] = out[6] = out[7] = static_cast<U>(0);

    return out;
  }

  // copy "const tensor2d" -> "tensor2s" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2s<U> () const
  {
    tensor2s<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[3] = static_cast<U>( m_data[1] );
    out[5] = static_cast<U>( m_data[2] );

    out[1] = out[2] = out[4] = static_cast<U>(0);

    return out;
  }

  // copy constructors : copy to std::vector
  // ---------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> () const
  {
    std::vector<U> out(3);

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // dimensions
  // ----------

  size_t size() const { return 3; }
  size_t ndim() const { return 3; }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];     }
  auto     begin() const { return &m_data[0];     }
  auto     end  () const { return &m_data[0] + 3; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

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

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { m_data[0] = m_data[1] = m_data[2] = D;                 }
  void setZero    (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); }
  void setOnes    (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); }
  void zeros      (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); }
  void ones       (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); }

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
  X           inline det   (                    ) const; // determinant
  tensor2d<X> inline inv   (                    ) const; // inverse

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; }
  const X& operator[](size_t i) const { return m_data[i]; }

  X&       operator()(size_t i, size_t j)
  {
    if (i == j) return m_data[i];
    else        return m_zero[0];
  }

  const X& operator()(size_t i, size_t j) const
  {
    if (i == j) return m_data[i];
    else        return m_zero[0];
  }

  // arithmetic operators: tensor2d ?= tensor2d
  // ------------------------------------------

  tensor2d<X>& operator*= (const tensor2d<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];

    return *this;
  }

  tensor2d<X>& operator+= (const tensor2d<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];

    return *this;
  }

  tensor2d<X>& operator-= (const tensor2d<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];

    return *this;
  }

  // arithmetic operators: tensor2d ?= tensor2
  // -----------------------------------------

  tensor2d<X>& operator*= (const tensor2 <X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[4];
    m_data[2] *= B[8];

    return *this;
  }

  tensor2d<X>& operator/= (const tensor2 <X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[4];
    m_data[2] /= B[8];

    return *this;
  }

  // arithmetic operators: tensor2d ?= tensor2s
  // ------------------------------------------

  tensor2d<X>& operator*= (const tensor2s<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[3];
    m_data[2] *= B[5];

    return *this;
  }

  tensor2d<X>& operator/= (const tensor2s<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[3];
    m_data[2] /= B[5];

    return *this;
  }

  // arithmetic operators: tensor2d ?= scalar
  // ----------------------------------------

  tensor2d<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;

    return *this;
  }

  tensor2d<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;

    return *this;
  }

  tensor2d<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;

    return *this;
  }

  tensor2d<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;

    return *this;
  }

  // equality operators: tensor2d == tensor2?
  // ----------------------------------------

  bool operator== ( const tensor2d<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = 0 ; j < 3 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  }

  bool operator== ( const tensor2s<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i ) {
      for ( size_t j = i ; j < 3 ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  }

}; // class tensor2d

// arithmetic operators: tensor2d = tensor2d ? tensor2d
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

template <class X> tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

template <class X> tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// arithmetic operators: tensor2d = tensor2d ? tensor2
// ---------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[4];
  C[2] = A[2] * B[8];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[4];
  C[2] = A[2] / B[8];

  return C;
}

// arithmetic operators: tensor2d = tensor2d ? tensor2s
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];
  C[2] = A[2] * B[5];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];
  C[2] = A[2] / B[5];

  return C;
}

// arithmetic operators: tensor2d = tensor2d ? scalar
// --------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// arithmetic operators: tensor2d = tensor2 ? tensor2d
// ---------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[4] * B[1];
  C[2] = A[8] * B[2];

  return C;
}

// arithmetic operators: tensor2d = tensor2s ? tensor2d
// ----------------------------------------------------

template <class X> tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];
  C[2] = A[5] * B[2];

  return C;
}


// arithmetic operators: tensor2d = scalar ? tensor2d
// --------------------------------------------------

template <class X> tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// =================================================================================================
// cppmat::cartesian3d::vector
// =================================================================================================

template<class X> class vector
{
private:

  // local data storage
  // - local container with the data
  X m_container[3];
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data = &m_container[0];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is initialized
  vector(){}

  // allocate tensor, initialize to constant "D"
  vector(X D) { m_data[0] = m_data[1] = m_data[2] = D; }

  // copy from raw pointer, correct storage
  vector(const X *D)
  {
    m_data[0] = D[0];
    m_data[1] = D[1];
    m_data[2] = D[2];
  }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(X *D) { m_data = D; }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(const X *D)
  {
    m_data[0] = D[0];
    m_data[1] = D[1];
    m_data[2] = D[2];
  }

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "vector" -> "vector" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator vector<U> () const
  {
    vector<U> out;

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // copy constructors : copy to std::vector
  // ---------------------------------------

  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator std::vector<U> () const
  {
    std::vector<U> out(3);

    out[0] = static_cast<U>( m_data[0] );
    out[1] = static_cast<U>( m_data[1] );
    out[2] = static_cast<U>( m_data[2] );

    return out;
  }

  // copy constructors : copy to Eigen object
  // ----------------------------------------

  // column
  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,3,1,Eigen::ColMajor> () const
  {
    Eigen::Matrix<U,3,1,Eigen::ColMajor> out;

    out(0) = static_cast<U>( m_data[0] );
    out(1) = static_cast<U>( m_data[1] );
    out(2) = static_cast<U>( m_data[2] );

    return out;
  }
  #endif

  // 'list'
  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,1,3,Eigen::RowMajor> () const
  {
    Eigen::Matrix<U,1,3,Eigen::RowMajor> out;

    out(0) = static_cast<U>( m_data[0] );
    out(1) = static_cast<U>( m_data[1] );
    out(2) = static_cast<U>( m_data[2] );

    return out;
  }
  #endif

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return 3; }
  size_t ndim() const { return 3; }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(1,3);

    return out;
  }

  // storage strides, optionally in bytes
  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> out = { 1 };

    out[0] *= sizeof(X);

    return out;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];     }
  auto     begin() const { return &m_data[0];     }
  auto     end  () const { return &m_data[0] + 3; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

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

  // length
  // ------

  X length() const
  {
    X C;

    C  = std::pow(m_data[0],2.);
    C += std::pow(m_data[1],2.);
    C += std::pow(m_data[2],2.);

    return std::sqrt(C);
  }

  // normalize to unit length
  // ------------------------

  void setUnitLength()
  {
    X C = length();

    if ( C <= static_cast<X>(0) ) return;

    m_data[0] /= C;
    m_data[1] /= C;
    m_data[2] /= C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { m_data[0] = m_data[1] = m_data[2] = D;                 }
  void setZero    (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); }
  void setOnes    (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); }
  void zeros      (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0); }
  void ones       (   ) { m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1); }

  // tensor products / operations
  // ----------------------------

  X          inline dot   (const vector  <X> &B) const; // dot    product: C   = A_i*B_i
  vector <X> inline dot   (const tensor2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2<X> inline dyadic(const vector  <X> &B) const; // dyadic product: C_ij = A_i*B_j
  vector <X> inline cross (const vector  <X> &B) const; // cross  product

  // index operators
  // ---------------

  X&       operator[](size_t i)       { return m_data[i]; }
  const X& operator[](size_t i) const { return m_data[i]; }
  X&       operator()(size_t i)       { return m_data[i]; }
  const X& operator()(size_t i) const { return m_data[i]; }

  // arithmetic operators: vector ?= vector
  // --------------------------------------

  vector<X>& operator*= (const vector<X> &B)
  {
    m_data[0] *= B[0];
    m_data[1] *= B[1];
    m_data[2] *= B[2];

    return *this;
  }

  vector<X>& operator/= (const vector<X> &B)
  {
    m_data[0] /= B[0];
    m_data[1] /= B[1];
    m_data[2] /= B[2];

    return *this;
  }

  vector<X>& operator+= (const vector<X> &B)
  {
    m_data[0] += B[0];
    m_data[1] += B[1];
    m_data[2] += B[2];

    return *this;
  }

  vector<X>& operator-= (const vector<X> &B)
  {
    m_data[0] -= B[0];
    m_data[1] -= B[1];
    m_data[2] -= B[2];

    return *this;
  }

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<X>& operator*= (const X &B)
  {
    m_data[0] *= B;
    m_data[1] *= B;
    m_data[2] *= B;

    return *this;
  }

  vector<X>& operator/= (const X &B)
  {
    m_data[0] /= B;
    m_data[1] /= B;
    m_data[2] /= B;

    return *this;
  }

  vector<X>& operator+= (const X &B)
  {
    m_data[0] += B;
    m_data[1] += B;
    m_data[2] += B;

    return *this;
  }

  vector<X>& operator-= (const X &B)
  {
    m_data[0] -= B;
    m_data[1] -= B;
    m_data[2] -= B;

    return *this;
  }

  // equality operator: vector == vector
  // -----------------------------------

  bool operator== ( const vector<X> &B )
  {
    for ( size_t i = 0 ; i < 3 ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

}; // class vector

// arithmetic operators: vector = vector ? vector
// ----------------------------------------------

template <class X> vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

template <class X> vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

template <class X> vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

template <class X> vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// arithmetic operators: vector = vector ? scalar
// ----------------------------------------------

template <class X> vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

template <class X> vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

template <class X> vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

template <class X> vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// arithmetic operators: vector = scalar ? vector
// ----------------------------------------------

template <class X> vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

template <class X> vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

template <class X> vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

template <class X> vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// =================================================================================================
// print to screen
// =================================================================================================

template<class X> void inline tensor4<X>::printf(std::string fmt) const
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
std::ostream& operator<<(std::ostream& out, tensor4<X>& src)
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          out << "(" << i << "," << j << "," << k << "," << l << ") " << src(i,j,k,l) << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[3],m_data[4],m_data[5]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[6],m_data[7],m_data[8]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2s<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[1],m_data[3],m_data[4]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[2],m_data[4],m_data[5]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2s<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline tensor2d<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_zero[0],m_zero[0]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_zero[0],m_data[1],m_zero[0]);
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_zero[0],m_zero[0],m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, tensor2d<X>& src)
{
  out << src(0,0) << ", " << src(0,1) << ", " << src(0,2) << ";" << std::endl;
  out << src(1,0) << ", " << src(1,1) << ", " << src(1,2) << ";" << std::endl;
  out << src(2,0) << ", " << src(2,1) << ", " << src(2,2) << ";" << std::endl;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X> void inline vector<X>::printf(std::string fmt) const
{
  std::printf((fmt+","+fmt+","+fmt+";\n").c_str(),m_data[0],m_data[1],m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template <class X>
std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  out << src(0) << ", " << src(1) << ", " << src(2) << ";" << std::endl;

  return out;
}

// =================================================================================================
// identity tensors
// =================================================================================================

tensor4<double> inline identity4()
{
  tensor4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==l and j==k )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4rt()
{
  tensor4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==k and j==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identityII()
{
  tensor4<double> I(0.0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          if ( i==j and k==l )
            I(i,j,k,l) = 1.;

  return I;
}

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4s()
{ return (identity4()+identity4rt())/2.; }

// -------------------------------------------------------------------------------------------------

tensor4<double> inline identity4d()
{ return identity4s()-identityII()/static_cast<double>(3); }

// -------------------------------------------------------------------------------------------------

tensor2d<double> inline identity2()
{
  tensor2d<double> I(1.);

  return I;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X> tensor4<X> inline tensor4<X>::ddot(const tensor4<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

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

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2s<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j) += (*this)(i,j,k,l)*B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor4<X>::ddot(const tensor2d<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j) += (*this)(i,j,k,k)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[4]*B[1] + m_data[8]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(k,l) += (*this)(i,j)*B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C += (*this)(i,j)*B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::ddot(const tensor2s<X> &B) const
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

template<class X> X inline tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[3]*B[1] + m_data[5]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[4] + m_data[2]*B[8];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[3] + m_data[2]*B[5];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[1] + m_data[2]*B[2];
}


// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::dot(const tensor2d<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2<X>::dot(const vector<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,k) += (*this)(i,j)*B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2s<X>::dot(const tensor2d<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i,j)*B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2s<X>::dot(const vector<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i) += (*this)(i,j)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2d<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,k) += (*this)(i,i)*B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::dot(const tensor2d<X> &B) const
{
  tensor2d<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i,i) += (*this)(i,i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline tensor2d<X>::dot(const vector<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i) += (*this)(i,i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2s<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(j) += (*this)(i)*B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::dot(const tensor2d<X> &B) const
{
  vector<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    C(i) += (*this)(i)*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline vector<X>::dot(const vector<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i=0; i<3; ++i )
    C += (*this)(i)*B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,k,l) += (*this)(i,j)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        C(i,j,k,k) += (*this)(i,j)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      for ( size_t l=0; l<3; ++l )
        C(i,i,k,l) += (*this)(i,i)*B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t k=0; k<3; ++k )
      C(i,i,k,k) += (*this)(i,i)*B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline vector<X>::dyadic(const vector<X> &B) const
{
  tensor2<X> C(static_cast<X>(0));

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      C(i,j) += (*this)(i)*B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> vector<X> inline vector<X>::cross(const vector<X> &B) const
{
  vector<X> C;

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
  tensor4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::RT() const
{
  tensor4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor4<X> inline tensor4<X>::LT() const
{
  tensor4<X> C;

  for ( size_t i=0; i<3; ++i )
    for ( size_t j=0; j<3; ++j )
      for ( size_t k=0; k<3; ++k )
        for ( size_t l=0; l<3; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::T() const
{
  tensor2<X> C;

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

template<class X> tensor2s<X> inline tensor2s<X>::T() const
{
  tensor2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];
  C[3] = m_data[3];
  C[4] = m_data[4];
  C[5] = m_data[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::T() const
{
  tensor2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// =================================================================================================
// miscellaneous
// =================================================================================================

template<class X> X inline tensor2<X>::trace() const
{
  return m_data[0]+m_data[4]+m_data[8];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::trace() const
{
  return m_data[0]+m_data[3]+m_data[5];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::trace() const
{
  return m_data[0]+m_data[1]+m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2<X>::det() const
{
  return ( m_data[0] * m_data[4] * m_data[8] +
           m_data[1] * m_data[5] * m_data[6] +
           m_data[2] * m_data[3] * m_data[7] ) -
         ( m_data[2] * m_data[4] * m_data[6] +
           m_data[1] * m_data[3] * m_data[8] +
           m_data[0] * m_data[5] * m_data[7] );
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2s<X>::det() const
{
  return (                     m_data[0] * m_data[3] * m_data[5] +
           static_cast<X>(2) * m_data[1] * m_data[2] * m_data[4] ) -
         (                     m_data[4] * m_data[4] * m_data[0] +
                               m_data[2] * m_data[2] * m_data[3] +
                               m_data[1] * m_data[1] * m_data[5] );
}

// -------------------------------------------------------------------------------------------------

template<class X> X inline tensor2d<X>::det() const
{
  return m_data[0] * m_data[1] * m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2<X> inline tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2<X> C;

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

template<class X> tensor2s<X> inline tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2s<X> C;

  C[0] = (m_data[3]*m_data[5]-m_data[4]*m_data[4]) / D;
  C[1] = (m_data[2]*m_data[4]-m_data[1]*m_data[5]) / D;
  C[2] = (m_data[1]*m_data[4]-m_data[2]*m_data[3]) / D;
  C[3] = (m_data[0]*m_data[5]-m_data[2]*m_data[2]) / D;
  C[4] = (m_data[2]*m_data[1]-m_data[0]*m_data[4]) / D;
  C[5] = (m_data[0]*m_data[3]-m_data[1]*m_data[1]) / D;
  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X> tensor2d<X> inline tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C;

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

} // namespace cartesian3d
} // namespace cppmat

#endif

