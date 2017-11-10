/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR_H
#define CPPMAT_TENSOR_H

#include "macros.h"

namespace cppmat {
namespace cartesian {

// =================================================================================================
// forward declaration
// =================================================================================================

template<class X> class tensor4;
template<class X> class tensor2;
template<class X> class tensor2s;
template<class X> class tensor2d;
template<class X> class vector;

// =================================================================================================
// cppmat::cartesian::tensor4
// =================================================================================================

template<class X> class tensor4
{
private:

  // local data storage
  // - local container with the data
  std::vector<X> m_container;
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data;
  // - number of dimensions
  size_t m_nd=0;
  // - size of the data array
  size_t m_size=0;

public:

  // constructors
  // ------------

  // allocate tensor, nothing is allocated
  tensor4(){}

  // allocate tensor, nothing is initialized
  tensor4(size_t nd) { resize(nd); }

  // allocate tensor, initialize to constant "D"
  tensor4(size_t nd, X D) { resize(nd); for ( size_t i=0; i<m_size; ++i ) m_data[i]=D; }

  // copy from raw pointer
  tensor4(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t nd, X *D)
  {
    // - change size settings
    m_nd   = nd;
    m_size = m_nd*m_nd*m_nd*m_nd;

    // - point to input pointer
    m_data = D;
  }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // resize container
  // ----------------

  void resize(size_t nd)
  {
    m_nd   = nd;
    m_size = m_nd*m_nd*m_nd*m_nd;

    m_container.resize(m_size);

    m_data = &m_container[0];
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
    tensor4<U> out(m_nd);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

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
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return m_size; }
  size_t ndim() const { return m_nd;   }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(4,m_nd);

    return out;
  }

  // storage strides, optionally in bytes
  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> out = { m_nd*m_nd*m_nd, m_nd*m_nd, m_nd, 1 };

    out[0] *= sizeof(X);
    out[1] *= sizeof(X);
    out[2] *= sizeof(X);
    out[3] *= sizeof(X);

    return out;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::abs(m_data[i]);

    return C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

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
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; }

  const X& operator()(size_t i, size_t j, size_t k, size_t l) const
  { return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l]; }

  // arithmetic operators: tensor4 ?= tensor4
  // ----------------------------------------

  tensor4<X>& operator*= (const tensor4<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i=0; i<m_size; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  tensor4<X>& operator/= (const tensor4<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i=0; i<m_size; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  tensor4<X>& operator+= (const tensor4<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i=0; i<m_size; ++i )
      m_data[i] += B[i];

    return *this;
  }

  tensor4<X>& operator-= (const tensor4<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i=0; i<m_size; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor4 ?= scalar
  // ---------------------------------------

  tensor4<X>& operator*= (const X &B)
  {
    for ( size_t i=0; i<m_size; ++i )
      m_data[i] *= B;

    return *this;
  }

  tensor4<X>& operator/= (const X &B)
  {
    for ( size_t i=0; i<m_size; ++i )
      m_data[i] /= B;

    return *this;
  }

  tensor4<X>& operator+= (const X &B)
  {
    for ( size_t i=0; i<m_size; ++i )
      m_data[i] += B;

    return *this;
  }

  tensor4<X>& operator-= (const X &B)
  {
    for ( size_t i=0; i<m_size; ++i )
      m_data[i] -= B;

    return *this;
  }

  // equality operators: tensor4 == tensor4
  // --------------------------------------

  bool operator== ( const tensor4<X> &B )
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i<size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

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
// cppmat::cartesian::tensor2
// =================================================================================================

template<class X> class tensor2
{
private:

  // local data storage
  // - local container with the data
  std::vector<X> m_container;
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data;
  // - number of dimensions
  size_t m_nd=0;
  // - size of the data array
  size_t m_size=0;

public:

  // constructors
  // ------------

  // allocate tensor, nothing is allocated
  tensor2(){}

  // allocate tensor, nothing is initialized
  tensor2(size_t nd) { resize(nd); }

  // allocate tensor, initialize to constant "D"
  tensor2(size_t nd, X D) { resize(nd); for ( size_t i=0; i<m_size; ++i ) m_data[i]=D; }

  // copy from raw pointer, correct storage
  tensor2(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  tensor2(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  tensor2(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        m_data[i*m_nd+j] = D(i,j);
  }
  #endif

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t nd, X *D)
  {
    // - change size settings
    m_nd   = nd;
    m_size = m_nd*m_nd;

    // - point to input pointer
    m_data = D;
  }

  // pointer to Eigen array
  // N.B. only possible for matching 'RowMajor' storage, the user has to keep the pointer alive
  #ifdef CPPMAT_EIGEN
  void map(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - update size settings
    m_nd   = static_cast<size_t>(D.rows());
    m_size = m_nd*m_nd;

    // - point to input array
    m_data = const_cast<X*>(&D(0));
  }
  #endif

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        m_data[i*m_nd+j] = D(i,j);
  }
  #endif

  // resize container
  // ----------------

  void resize(size_t nd)
  {
    m_nd   = nd;
    m_size = m_nd*m_nd;

    m_container.resize(m_size);

    m_data = &m_container[0];
  }

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

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
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy constructors : copy to Eigen object
  // ----------------------------------------

  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> () const
  {
    Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> out(m_nd,m_nd);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out(i) = static_cast<U>( m_data[i] );

    return out;
  }
  #endif

  // conversions : convert to incompatible tensor objects -> SOME TERMS ARE DISCARDED
  // --------------------------------------------------------------------------------

  // "tensor2 -> tensor2s" (output is symmetrized: "out(i,j) = ( this(i,j) + this(j,i) ) / 2.")
  tensor2s<X> astensor2s()
  {
    tensor2s<X> out(m_nd);

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i ; j < m_nd ; ++j )
        out[i*m_nd-(i-1)*i/2+j-i] = ( m_data[i*m_nd+j] + m_data[j*m_nd+i] ) / static_cast<X>(2);

    return out;
  }

  // "tensor2 -> tensor2d" (all off-diagonal terms are discarded)
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out(m_nd);

    for ( size_t i = 0 ; i < m_nd ; ++i )
      out[i] = m_data[i*m_nd+i];

    return out;
  }

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return m_size; }
  size_t ndim() const { return m_nd;   }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(2,m_nd);

    return out;
  }

  // storage strides, optionally in bytes
  std::vector<size_t> strides(bool bytes=false) const
  {
    std::vector<size_t> out = { m_nd, 1 };

    out[0] *= sizeof(X);
    out[1] *= sizeof(X);

    return out;
  }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::abs(m_data[i]);

    return C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

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

  X&       operator[](size_t i          )       { return m_data[i];        }
  const X& operator[](size_t i          ) const { return m_data[i];        }
  X&       operator()(size_t i, size_t j)       { return m_data[i*m_nd+j]; }
  const X& operator()(size_t i, size_t j) const { return m_data[i*m_nd+j]; }

  // arithmetic operators: tensor2 ?= tensor2
  // ----------------------------------------

  tensor2<X>& operator*= (const tensor2<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  tensor2<X>& operator/= (const tensor2<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  tensor2<X>& operator+= (const tensor2<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  tensor2<X>& operator-= (const tensor2<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor2 ?= tensor2s
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        // - extract value
        X b = B[i*m_nd-(i-1)*i/2+j-i];
        // - store symmetrically
                      m_data[i*m_nd+j] *= b;
        if ( i != j ) m_data[j*m_nd+i] *= b;
      }
    }

    return *this;
  }

  tensor2<X>& operator/= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        // - extract value
        X b = B[i*m_nd-(i-1)*i/2+j-i];
        // - store symmetrically
                      m_data[i*m_nd+j] /= b;
        if ( i != j ) m_data[j*m_nd+i] /= b;
      }
    }

    return *this;
  }

  tensor2<X>& operator+= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        // - extract value
        X b = B[i*m_nd-(i-1)*i/2+j-i];
        // - store symmetrically
                      m_data[i*m_nd+j] += b;
        if ( i != j ) m_data[j*m_nd+i] += b;
      }
    }

    return *this;
  }

  tensor2<X>& operator-= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        // - extract value
        X b = B[i*m_nd-(i-1)*i/2+j-i];
        // - store symmetrically
                      m_data[i*m_nd+j] -= b;
        if ( i != j ) m_data[j*m_nd+i] -= b;
      }
    }

    return *this;
  }

  // arithmetic operators: tensor2 ?= tensor2d
  // -----------------------------------------

  tensor2<X>& operator*= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = 0 ; j < m_nd ; ++j ) {
        if ( i == j ) m_data[i*m_nd+i] *= B[i];
        else          m_data[i*m_nd+j]  = static_cast<X>(0);
      }
    }

    return *this;
  }

  tensor2<X>& operator+= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i*m_nd+i] += B[i];

    return *this;
  }

  tensor2<X>& operator-= (const tensor2d<X> &B)
  {
    assert( ndim() == B.ndim() );

    for ( size_t i=0; i<m_nd; ++i )
      m_data[i*m_nd+i] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor2 ?= scalar
  // ---------------------------------------

  tensor2<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  tensor2<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  tensor2<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  tensor2<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B;

    return *this;
  }

  // equality operators: tensor2 == tensor2?
  // ---------------------------------------

  bool operator== ( const tensor2<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0;  i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2s<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( m_data[i*m_nd+j] != B(i,j) )
          return false;

    return true;
  }

  bool operator== ( const tensor2d<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( m_data[i*m_nd+j] != B(i,j) )
          return false;

    return true;
  }

  // structure-check operators
  // -------------------------

  bool issymmetric()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i+1 ; j < m_nd ; ++j )
        if ( m_data[ i*m_nd + j ] != m_data[ j*m_nd + i ] )
          return false;

    return true;
  }

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = 0 ; j < m_nd ; ++j )
        if ( i != j )
          if ( m_data[ i*m_nd + j ] )
            return false;

    return true;
  }

}; // class tensor2

// arithmetic operators: tensor2 = tensor2 ? tensor2
// -------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// arithmetic operators: tensor2 = tensor2 ? tensor2s
// --------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i*nd + j ] + B[ i ];
      else          C[ i*nd + j ] = A[ i*nd + j ];
    }
  }

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

template <class X> tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

template <class X> tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

template <class X> tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// arithmetic operators: tensor2 = tensor2s ? tensor2
// --------------------------------------------------

template <class X> tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

template <class X> tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i ] + B[ i*nd + j ];
      else          C[ i*nd + j ] =          B[ i*nd + j ];
    }
  }

  return C;
}

template <class X> tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

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

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

template <class X> tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

template <class X> tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

template <class X> tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// cppmat::cartesian::tensor2s (symmetric storage of "cppmat::tensor")
// =================================================================================================

template<class X> class tensor2s
{
private:

  // local data storage
  // - local container with the data
  std::vector<X> m_container;
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data;
  // - number of dimensions
  size_t m_nd=0;
  // - size of the data array
  size_t m_size=0;

public:

  // constructors
  // ------------

  // allocate tensor, nothing is allocated
  tensor2s(){}

  // allocate tensor, nothing is initialized
  tensor2s(size_t nd) { resize(nd); }

  // allocate tensor, initialize to constant "D"
  tensor2s(size_t nd, X D) { resize(nd); for ( size_t i=0; i<m_size; ++i ) m_data[i]=D; }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t nd, X *D)
  {
    m_nd   = nd;
    m_size = (m_nd+1)*m_nd/2;

    m_data = D;
  }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // raw pointer, stored as if it was "tensor2"
  void copyDense(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = i+1 ; j < m_nd ; ++j )
          assert( D[i*m_nd+j] == D[j*m_nd+i] );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < nd ; ++i )
      for ( size_t j = i ; j < nd ; ++j )
        m_data[i*m_nd-(i-1)*i/2+j-i] = D[i*m_nd+j];
  }

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = i+1 ; j < m_nd ; ++j )
          assert( D(i,j) == D(j,i) );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i ; j < m_nd ; ++j )
        m_data[i*m_nd-(i-1)*i/m_nd+j-i] = D(i,j);
  }
  #endif

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - check the input for symmetry (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = i+1 ; j < m_nd ; ++j )
          assert( D(i,j) == D(j,i) );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i ; j < m_nd ; ++j )
        m_data[i*m_nd-(i-1)*i/m_nd+j-i] = D(i,j);
  }
  #endif

  // resize container
  // ----------------

  void resize(size_t nd)
  {
    m_nd   = nd;
    m_size = (m_nd+1)*m_nd/2;

    m_container.resize(m_size);

    m_data = &m_container[0];
  }

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2s" -> "tensor2s" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2s<U> () const
  {
    tensor2s<U> out(m_nd);

    for ( size_t i=0; i<size(); ++i )
      out[i] = static_cast<U>(m_data[i]);

    return out;
  }

  // copy "const tensor2s" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        // - get item
        U b = static_cast<U>(m_data[i*m_nd-(i-1)*i/2+j-i]);
        // - store item, and symmetric copy
        out[i*m_nd+j] = b;
        out[j*m_nd+i] = b;
      }
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
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // conversions : convert to incompatible tensor objects -> SOME TERMS ARE DISCARDED
  // --------------------------------------------------------------------------------

  // "tensor2s -> tensor2d" (all off-diagonal terms are discarded)
  tensor2d<X> astensor2d()
  {
    tensor2d<X> out(m_nd);

    for ( size_t i=0; i<m_nd; ++i )
      out[i] = m_data[ i*m_nd - (i-1)*i/2 ];

    return out;
  }

  // dimensions
  // ----------

  size_t size() const { return m_size; }
  size_t ndim() const { return m_nd;   }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::abs(m_data[i]);

    return C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

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

  X&       operator[](size_t i )       { return m_data[i]; }
  const X& operator[](size_t i ) const { return m_data[i]; }

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
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  tensor2s<X>& operator/= (const tensor2s<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  tensor2s<X>& operator+= (const tensor2s<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  tensor2s<X>& operator-= (const tensor2s<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor2s ?= tensor2d
  // ------------------------------------------

  tensor2s<X>& operator*= (const tensor2d<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j=i; j<m_nd; ++j ) {
        if ( i == j ) m_data[i*m_nd-(i-1)*i/2    ] *= B[i];
        else          m_data[i*m_nd-(i-1)*i/2+j-i]  = static_cast<X>(0);
      }
    }

    return *this;
  }

  tensor2s<X>& operator+= (const tensor2d<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i*m_nd-(i-1)*i/2] += B[i];

    return *this;
  }

  tensor2s<X>& operator-= (const tensor2d<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i*m_nd-(i-1)*i/2] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor2s ?= scalar
  // ----------------------------------------

  tensor2s<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  tensor2s<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  tensor2s<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  tensor2s<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B;

    return *this;
  }

  // equality operators: tensor2 == tensor2?
  // ---------------------------------------

  bool operator== ( const tensor2s<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0; i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
        if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(i,j) ) return false;
        if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(j,i) ) return false;
      }
    }

    return true;
  }

  bool operator== ( const tensor2d<X> &B )
  {
    assert( ndim() == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i ; j < m_nd ; ++j )
        if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(i,j) ) return false;

    return true;
  }

  // structure-check operators
  // -------------------------

  bool isdiagonal()
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      for ( size_t j = i+1 ; j < m_nd ; ++j )
        if ( m_data[i*m_nd-(i-1)*i/2+j-i] )
          return false;

    return true;
  }

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

  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
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
  size_t nd = A.ndim();
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
  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
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
// cppmat::cartesian::tensor2d (symmetric storage of "cppmat::tensor::tensor")
// =================================================================================================

template<class X> class tensor2d
{
private:

  // local data storage
  // - local container with the data
  std::vector<X> m_container;
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data;
  // - number of dimensions
  size_t m_nd=0;
  // - size of the data array
  size_t m_size=0;

  // dummy parameter, used to return "0" for any off-diagonal entry
  X m_zero[1];

public:

  // constructors
  // ------------

  // allocate tensor, nothing is allocated (but the dummy parameter is set)
  tensor2d(){ m_zero[0] = static_cast<X>(0); }

  // allocate tensor, nothing is initialized (but the dummy parameter is set)
  tensor2d(size_t nd) { m_zero[0] = static_cast<X>(0); resize(nd); }

  // allocate tensor, initialize to constant "D" (and set the dummy parameter)
  tensor2d(size_t nd, X D)
  {
    m_zero[0] = static_cast<X>(0);

    resize(nd);

    for ( size_t i=0; i<m_size; ++i ) m_data[i]=D;
  }

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t nd, X *D)
  {
    m_nd   = nd;
    m_size = m_nd;

    m_data = D;
  }

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // raw pointer, stored as if it was "tensor2"
  void copyDense(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = 0 ; j < m_nd ; ++j )
          if ( i != j )
            assert( ! D[i*m_nd+j] );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < nd ; ++i )
      m_data[i] = D[i*m_nd+i];
  }

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = 0 ; j < m_nd ; ++j )
          if ( i != j )
            assert( ! D(i,j) );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] = D(i,i);
  }
  #endif

  // Eigen array, stored as if it was "tensor2"
  #ifdef CPPMAT_EIGEN
  void copyDense(const Eigen::Matrix<X,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &D)
  {
    // - check size
    assert( D.rows() == D.cols() );

    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - check that input is diagonal (code eliminated if "NDEBUG" is defined)
    #ifndef NDEBUG
      for ( size_t i = 0 ; i < m_nd ; ++i )
        for ( size_t j = 0 ; j < m_nd ; ++j )
          if ( i != j )
            assert( ! D(i,j) );
    #endif

    // - copy diagonal
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] = D(i,i);
  }
  #endif

  // resize container
  // ----------------

  void resize(size_t nd)
  {
    m_nd   = nd;
    m_size = m_nd;

    m_container.resize(m_size);

    m_data = &m_container[0];
  }

  // copy constructors : copy to a new tensor object
  // -----------------------------------------------

  // copy "tensor2d" -> "tensor2d" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2d<U> () const
  {
    tensor2d<U> out(m_nd);

    for ( size_t i = 0 ; i < m_nd ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy "const tensor2d" -> "tensor2" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2<U> () const
  {
    tensor2<U> out(m_nd,static_cast<U>(0));

    for ( size_t i = 0 ; i < m_nd ; ++i )
      out[i*m_nd+i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy "const tensor2d" -> "tensor2s" ( + change of type )
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator tensor2s<U> () const
  {
    tensor2s<U> out(m_nd,static_cast<U>(0));

    for ( size_t i = 0 ; i < m_nd ; ++i )
      out[i*m_nd-(i-1)*i/2] = static_cast<U>( m_data[i] );

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
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // dimensions
  // ----------

  size_t size() const { return m_size; }
  size_t ndim() const { return m_nd;   }

  // pointer / iterators
  // -------------------

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::abs(m_data[i]);

    return C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

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
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  tensor2d<X>& operator+= (const tensor2d<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  tensor2d<X>& operator-= (const tensor2d<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: tensor2d ?= tensor2
  // -----------------------------------------

  tensor2d<X>& operator*= (const tensor2 <X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd; ++i )
      m_data[i] *= B[i*m_nd+i];

    return *this;
  }

  tensor2d<X>& operator/= (const tensor2 <X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd; ++i )
      m_data[i] /= B[i*m_nd+i];

    return *this;
  }

  // arithmetic operators: tensor2d ?= tensor2s
  // ------------------------------------------

  tensor2d<X>& operator*= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] *= B[i*m_nd-(i-1)*i/2];

    return *this;
  }

  tensor2d<X>& operator/= (const tensor2s<X> &B)
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] /= B[i*m_nd-(i-1)*i/2];

    return *this;
  }

  // arithmetic operators: tensor2d ?= scalar
  // ----------------------------------------

  tensor2d<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] *= B;

    return *this;
  }

  tensor2d<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] /= B;

    return *this;
  }

  tensor2d<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] += B;

    return *this;
  }

  tensor2d<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < m_nd ; ++i )
      m_data[i] -= B;

    return *this;
  }

  // equality operators: tensor2d == tensor2?
  // ----------------------------------------

  bool operator== ( const tensor2d<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

  bool operator== ( const tensor2<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = 0 ; j < m_nd ; ++j ) {
        if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
        else          { if (              B(i,j) ) return false; }
      }
    }

    return true;
  }

  bool operator== ( const tensor2s<X> &B )
  {
    assert( m_nd == B.ndim() );

    for ( size_t i = 0 ; i < m_nd ; ++i ) {
      for ( size_t j = i ; j < m_nd ; ++j ) {
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

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] * B[ i*nd + i ];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i=0; i<nd; ++i )
    C[i] = A[i] * B[ i*nd - (i-1)*i/2 ];

  return C;
}

template <class X> tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
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

  size_t nd = A.ndim();
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
// cppmat::cartesian::vector
// =================================================================================================

template<class X> class vector
{
private:

  // local data storage
  // - local container with the data
  std::vector<X> m_container;
  // - pointer to the data container
  //   normally points to m_container, but may point to some external object using "map"
  X *m_data;
  // - number of dimensions
  size_t m_nd=0;
  // - size of the data array
  size_t m_size=0;

public:


  // constructors
  // ------------

  // allocate tensor, nothing is allocated
  vector(){}

  // allocate tensor, nothing is initialized
  vector(size_t nd) { resize(nd); }

  // allocate tensor, initialize to constant "D"
  vector(size_t nd, X D) { resize(nd); for ( size_t i=0; i<m_size; ++i ) m_data[i]=D; }

  // copy from raw pointer, correct storage
  vector(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  vector(const Eigen::Matrix<X,1,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - change size of internal storage
    resize(static_cast<size_t>(D.cols()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // copy from Eigen array
  #ifdef CPPMAT_EIGEN
  vector(const Eigen::Matrix<X,Eigen::Dynamic,1,Eigen::ColMajor> &D)
  {
    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // map external pointer
  // --------------------

  // raw pointer
  // N.B. the user is responsible for the correct storage and to keep the pointer alive
  void map(size_t nd, X *D)
  {
    // - change size settings
    m_nd   = nd;
    m_size = m_nd;

    // - point to input pointer
    m_data = D;
  }

  // pointer to Eigen array
  #ifdef CPPMAT_EIGEN
  void map(const Eigen::Matrix<X,1,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - update size settings
    m_nd   = static_cast<size_t>(D.cols());
    m_size = m_nd;

    // - point to input array
    m_data = const_cast<X*>(&D(0));
  }
  #endif

  // pointer to Eigen array
  #ifdef CPPMAT_EIGEN
  void map(const Eigen::Matrix<X,Eigen::Dynamic,1,Eigen::ColMajor> &D)
  {
    // - update size settings
    m_nd   = static_cast<size_t>(D.rows());
    m_size = m_nd;

    // - point to input array
    m_data = const_cast<X*>(&D(0));
  }
  #endif

  // copy from external data array
  // -----------------------------

  // raw pointer, correct storage
  void copy(size_t nd, const X *D)
  {
    // - change size of internal storage
    resize(nd);

    // - copy as is
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D[i];
  }

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,1,Eigen::Dynamic,Eigen::RowMajor> &D)
  {
    // - change size of internal storage
    resize(static_cast<size_t>(D.cols()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // Eigen array
  #ifdef CPPMAT_EIGEN
  void copy(const Eigen::Matrix<X,Eigen::Dynamic,1,Eigen::ColMajor> &D)
  {
    // - change size of internal storage
    resize(static_cast<size_t>(D.rows()));

    // - copy from input
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] = D(i);
  }
  #endif

  // resize container
  // ----------------

  void resize(size_t nd)
  {
    m_nd   = nd;
    m_size = m_nd;

    m_container.resize(m_size);

    m_data = &m_container[0];
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
    vector<U> out(m_nd);

    for ( size_t i = 0 ; i < size() ; ++i )
      out[i] = static_cast<U>( m_data[i] );

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
    std::vector<U> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out[i] = static_cast<U>( m_data[i] );

    return out;
  }

  // copy constructors : copy to Eigen object
  // ----------------------------------------

  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,1,Eigen::Dynamic,Eigen::RowMajor> () const
  {
    Eigen::Matrix<U,1,Eigen::Dynamic,Eigen::RowMajor> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out(i) = static_cast<U>( m_data[i] );

    return out;
  }
  #endif

  #ifdef CPPMAT_EIGEN
  template<\
    typename U,typename V=X,\
    typename=typename std::enable_if<std::is_convertible<X,U>::value>::type\
  >
  operator Eigen::Matrix<U,Eigen::Dynamic,1,Eigen::ColMajor> () const
  {
    Eigen::Matrix<U,Eigen::Dynamic,1,Eigen::ColMajor> out(m_size);

    for ( size_t i = 0 ; i < m_size ; ++i )
      out(i) = static_cast<U>( m_data[i] );

    return out;
  }
  #endif

  // dimensions
  // ----------

  // dimensions
  size_t size() const { return m_size; }
  size_t ndim() const { return m_nd;   }

  // shape
  std::vector<size_t> shape() const
  {
    std::vector<size_t> out(1,m_nd);

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

  const X* data () const { return &m_data[0];          }
  auto     begin() const { return &m_data[0];          }
  auto     end  () const { return &m_data[0] + m_size; }

  // print to screen
  // ---------------

  // formatted print (code below); NB also "operator<<" is defined below
  void printf(std::string fmt) const;

  // norm
  // ----

  X norm() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::abs(m_data[i]);

    return C;
  }

  // length
  // ------

  X length() const
  {
    X C = static_cast<X>(0);

    for ( size_t i = 0 ; i < m_size ; ++i )
      C += std::pow(m_data[i],2.);

    return std::sqrt(C);
  }

  // normalize to unit length
  // ------------------------

  void setUnitLength()
  {
    X C = length();

    if ( C <= static_cast<X>(0) ) return;

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= C;
  }

  // initialize all entries to zero/one/constant
  // -------------------------------------------

  void setConstant(X D) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = D;                 }
  void setZero    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void setOnes    (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }
  void zeros      (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(0); }
  void ones       (   ) { for ( size_t i=0; i<m_size; ++i ) m_data[i] = static_cast<X>(1); }

  // tensor products / operations
  // ----------------------------

  X          inline dot   (const vector  <X> &B) const; // dot    product: C   = A_i*B_i
  vector <X> inline dot   (const tensor2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2<X> inline dyadic(const vector  <X> &B) const; // dyadic product: C_ij = A_i*B_j
  vector <X> inline cross (const vector  <X> &B) const; // cross  product (only in 3D)

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
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B[i];

    return *this;
  }

  vector<X>& operator/= (const vector<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B[i];

    return *this;
  }

  vector<X>& operator+= (const vector<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B[i];

    return *this;
  }

  vector<X>& operator-= (const vector<X> &B)
  {
    assert( m_size == B.size() );
    assert( m_nd   == B.ndim() );

    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B[i];

    return *this;
  }

  // arithmetic operators: vector ?= scalar
  // --------------------------------------

  vector<X>& operator*= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] *= B;

    return *this;
  }

  vector<X>& operator/= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] /= B;

    return *this;
  }

  vector<X>& operator+= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] += B;

    return *this;
  }

  vector<X>& operator-= (const X &B)
  {
    for ( size_t i = 0 ; i < m_size ; ++i )
      m_data[i] -= B;

    return *this;
  }

  // equality operator: vector == vector
  // -----------------------------------

  bool operator== ( const vector<X> &B )
  {
    for ( size_t i = 0 ; i < size() ; ++i )
      if ( m_data[i] != B[i] )
        return false;

    return true;
  }

}; // class vector

// arithmetic operators: vector = vector ? vector
// ----------------------------------------------

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
// ----------------------------------------------

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
// ----------------------------------------------

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
      if (i == j) std::printf((fmt+",").c_str(),m_data[i]);
      else        std::printf((fmt+",").c_str(),m_zero[0]);
    }
    j = nd-1;
    if (i == j) std::printf((fmt+";\n").c_str(),m_data[i]);
    else        std::printf((fmt+";\n").c_str(),m_zero[0]);
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

} // namespace cartesian
} // namespace cppmat

#endif

