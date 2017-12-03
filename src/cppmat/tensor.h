/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR_H
#define CPPMAT_TENSOR_H

#include "macros.h"
#include "tensor2.h"
#include "tensor3.h"

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

template<class X>
class tensor4
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_nd=0;    // number of dimensions (rank == 4 -> shape == [m_nd,m_nd,m_nd,m_nd])
  size_t         m_size=0;  // == m_data.size() == m_nd*m_nd*m_nd*m_nd

public:

  // constructors
  tensor4(){};
  tensor4(size_t nd);
  tensor4(size_t nd, X D);

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j, size_t k, size_t l);
  const X& operator()(size_t i, size_t j, size_t k, size_t l) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  tensor4<X>& operator*= (const tensor4<X> &B); // A_ijkl *= B_ijkl
  tensor4<X>& operator/= (const tensor4<X> &B); // A_ijkl /= B_ijkl
  tensor4<X>& operator+= (const tensor4<X> &B); // A_ijkl += B_ijkl
  tensor4<X>& operator-= (const tensor4<X> &B); // A_ijkl -= B_ijkl
  tensor4<X>& operator*= (const         X  &B); // A_ijkl *= B
  tensor4<X>& operator/= (const         X  &B); // A_ijkl /= B
  tensor4<X>& operator+= (const         X  &B); // A_ijkl += B
  tensor4<X>& operator-= (const         X  &B); // A_ijkl -= B

  // tensor products / operations
  tensor4<X> inline ddot(const tensor4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2<X> inline ddot(const tensor2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> inline ddot(const tensor2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> inline ddot(const tensor2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor4<X> inline T   (                    ) const; // transposition   : C_lkji = A_ijkl
  tensor4<X> inline RT  (                    ) const; // transposition   : C_ijlk = A_ijkl
  tensor4<X> inline LT  (                    ) const; // transposition   : C_jikl = A_ijkl

  // equality operators
  bool operator== (const tensor4<X> &B) const; // A_ijkl == B_ijkl

  // basic algebra
  X norm() const; // sum(| A_ijkl |)

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian::tensor2
// =================================================================================================

template<class X>
class tensor2
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_nd=0;    // number of dimensions (rank == 2 -> shape == [m_nd,m_nd])
  size_t         m_size=0;  // == m_data.size() == m_nd*m_nd

public:

  // constructors
  tensor2(){};
  tensor2(size_t nd);
  tensor2(size_t nd, X D);

  // cast into another object
  template<class U> U cast() const;

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j);
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  tensor2<X>& operator*= (const tensor2 <X> &B); // A_ij *= B_ij
  tensor2<X>& operator/= (const tensor2 <X> &B); // A_ij /= B_ij
  tensor2<X>& operator+= (const tensor2 <X> &B); // A_ij += B_ij
  tensor2<X>& operator-= (const tensor2 <X> &B); // A_ij -= B_ij
  tensor2<X>& operator*= (const tensor2s<X> &B); // A_ij *= B_ij
  tensor2<X>& operator/= (const tensor2s<X> &B); // A_ij /= B_ij
  tensor2<X>& operator+= (const tensor2s<X> &B); // A_ij += B_ij
  tensor2<X>& operator-= (const tensor2s<X> &B); // A_ij -= B_ij
  tensor2<X>& operator*= (const tensor2d<X> &B); // A_ij *= B_ij
  tensor2<X>& operator+= (const tensor2d<X> &B); // A_ii += B_ii
  tensor2<X>& operator-= (const tensor2d<X> &B); // A_ii -= B_ii
  tensor2<X>& operator*= (const          X  &B); // A_ij *= B
  tensor2<X>& operator/= (const          X  &B); // A_ij /= B
  tensor2<X>& operator+= (const          X  &B); // A_ij += B
  tensor2<X>& operator-= (const          X  &B); // A_ij -= B

  // tensor products / operations
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
  tensor2<X> inline T     (                    ) const; // transpose       : C_ij   = A_ji
  X          inline trace (                    ) const; // trace           : A_ii
  X          inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // equality operators
  bool operator== (const tensor2 <X> &B) const; // A_ij == B_ij
  bool operator== (const tensor2s<X> &B) const; // A_ij == A_ij
  bool operator== (const tensor2d<X> &B) const; // A_ij == B_ij

  // structure check
  bool issymmetric() const; // A_ij == A_ji
  bool isdiagonal() const;  // A_ij == 0, i != j

  // basic algebra
  X norm() const; // sum(| A_ij |)

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian::tensor2s (symmetric storage of "cppmat::cartesian::tensor")
// =================================================================================================

template<class X>
class tensor2s
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_nd=0;    // number of dimensions (rank == 2 -> shape == [m_nd,m_nd])
  size_t         m_size=0;  // == m_data.size() == (m_nd+1)*m_nd/2

public:

  // constructors
  tensor2s(){};
  tensor2s(size_t nd);
  tensor2s(size_t nd, X D);

  // cast into another object
  template<class U> U cast() const;

  // automatic conversion to "tensor2"
  #ifndef CPPMAT_NOCONVERT
  operator tensor2<X> () const;
  #endif

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j);
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  tensor2s<X>& operator*= (const tensor2s<X> &B); // A_ij *= B_ij
  tensor2s<X>& operator/= (const tensor2s<X> &B); // A_ij /= B_ij
  tensor2s<X>& operator+= (const tensor2s<X> &B); // A_ij += B_ij
  tensor2s<X>& operator-= (const tensor2s<X> &B); // A_ij -= B_ij
  tensor2s<X>& operator*= (const tensor2d<X> &B); // A_ij *= B_ij
  tensor2s<X>& operator+= (const tensor2d<X> &B); // A_ii += B_ii
  tensor2s<X>& operator-= (const tensor2d<X> &B); // A_ii -= B_ii
  tensor2s<X>& operator*= (const          X  &B); // A_ij *= B
  tensor2s<X>& operator/= (const          X  &B); // A_ij /= B
  tensor2s<X>& operator+= (const          X  &B); // A_ij += B
  tensor2s<X>& operator-= (const          X  &B); // A_ij -= B

  // tensor products / operations
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
  tensor2s<X> inline T     (                    ) const; // transpose       : C_ij   = A_ji
  X           inline trace (                    ) const; // trace           : A_ii
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2s<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // equality operators
  bool operator== (const tensor2 <X> &B) const; // A_ij == A_ij
  bool operator== (const tensor2s<X> &B) const; // A_ij == B_ij
  bool operator== (const tensor2d<X> &B) const; // A_ij == B_ij

  // structure check
  bool isdiagonal() const; // A_ij == 0, i != j

  // basic algebra
  X norm() const; // sum(| A_ij |)

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian::tensor2d (symmetric storage of "cppmat::cartesian::tensor")
// =================================================================================================

template<class X>
class tensor2d
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_nd=0;    // number of dimensions (rank == 2 -> shape == [m_nd,m_nd])
  size_t         m_size=0;  // == m_data.size() == m_nd
  X              m_zero[1]; // dummy parameter, used to return "0" for any off-diagonal entry

public:

  // constructors
  tensor2d();
  tensor2d(size_t nd);
  tensor2d(size_t nd, X D);

  // cast into another object
  template<class U> U cast() const;

  // automatic conversion to "tensor2" or "tensor2s"
  #ifndef CPPMAT_NOCONVERT
  operator tensor2 <X> () const;
  operator tensor2s<X> () const;
  #endif

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j);
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  tensor2d<X>& operator*= (const tensor2d<X> &B); // A_ii * B_ii
  tensor2d<X>& operator+= (const tensor2d<X> &B); // A_ii + B_ii
  tensor2d<X>& operator-= (const tensor2d<X> &B); // A_ii - B_ii
  tensor2d<X>& operator*= (const tensor2 <X> &B); // A_ii * B_ii
  tensor2d<X>& operator/= (const tensor2 <X> &B); // A_ii / B_ii
  tensor2d<X>& operator*= (const tensor2s<X> &B); // A_ii * B_ii
  tensor2d<X>& operator/= (const tensor2s<X> &B); // A_ii / B_ii
  tensor2d<X>& operator*= (const          X  &B); // A_ii * B
  tensor2d<X>& operator/= (const          X  &B); // A_ii / B
  tensor2d<X>& operator+= (const          X  &B); // A_ii + B
  tensor2d<X>& operator-= (const          X  &B); // A_ii - B

  // tensor products / operations
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
  tensor2d<X> inline T     (                    ) const; // transpose      : C_ij   = A_ji
  X           inline trace (                    ) const; // trace          : A_ii
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2d<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // equality operators
  bool operator== (const tensor2 <X> &B) const; // A_ij == B_ij
  bool operator== (const tensor2s<X> &B) const; // A_ij == A_ij
  bool operator== (const tensor2d<X> &B) const; // A_ij == B_ij

  // basic algebra
  X norm() const; // sum(| A_ii |)

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian::vector
// =================================================================================================

template<class X>
class vector
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_nd=0;    // number of dimensions (rank == 1 -> shape == [m_nd])
  size_t         m_size=0;  // == m_data.size() == m_nd

public:

  // constructors
  vector(){};
  vector(size_t nd);
  vector(size_t nd, X D);

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i);
  const X& operator()(size_t i) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  vector<X>& operator*= (const vector<X> &B); // A_i * B_i
  vector<X>& operator/= (const vector<X> &B); // A_i / B_i
  vector<X>& operator+= (const vector<X> &B); // A_i + B_i
  vector<X>& operator-= (const vector<X> &B); // A_i - B_i
  vector<X>& operator*= (const        X  &B); // A_i * B
  vector<X>& operator/= (const        X  &B); // A_i / B
  vector<X>& operator+= (const        X  &B); // A_i + B
  vector<X>& operator-= (const        X  &B); // A_i - B

  // tensor products / operations
  X          inline dot   (const vector  <X> &B) const; // dot    product: C   = A_i*B_i
  vector <X> inline dot   (const tensor2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X> inline dot   (const tensor2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2<X> inline dyadic(const vector  <X> &B) const; // dyadic product: C_ij = A_i*B_j
  vector <X> inline cross (const vector  <X> &B) const; // cross  product (only in 3D)

  // equality operators
  bool operator== (const vector<X> &B) const; // A_i == B_i

  // basic algebra
  X norm() const;       // sum(| A_i |)
  X length() const;     // sqrt(sum(pow(A_i,2.)))
  void setUnitLength(); // A_i /= A.length()

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X> inline tensor4 <X> operator* (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator/ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator+ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator- (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator* (const tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator/ (const tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator+ (const tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator- (const tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator* (const          X  &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator/ (const          X  &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator+ (const          X  &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator- (const          X  &A, const tensor4 <X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator/ (const tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator* (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator* (const          X  &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const          X  &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const          X  &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const          X  &A, const tensor2 <X> &B);
template<class X> inline tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2s<X> operator* (const tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator/ (const tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator- (const tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator- (const tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator* (const          X  &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator/ (const          X  &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const          X  &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const          X  &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const          X  &A, const tensor2d<X> &B);
template<class X> inline tensor2s<X> operator- (const          X  &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2d<X> operator/ (const tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const          X  &A, const tensor2d<X> &B);
template<class X> inline vector  <X> operator* (const vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator/ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator+ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator- (const vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator* (const vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator/ (const vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator+ (const vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator- (const vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator* (const          X  &A, const vector  <X> &B);
template<class X> inline vector  <X> operator/ (const          X  &A, const vector  <X> &B);
template<class X> inline vector  <X> operator+ (const          X  &A, const vector  <X> &B);
template<class X> inline vector  <X> operator- (const          X  &A, const vector  <X> &B);

// =================================================================================================
// tensor products / operations
// =================================================================================================

template<class X> inline tensor4 <X> ddot  (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor2 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor2s<X> &A, const tensor4 <X> &B);
template<class X> inline tensor2 <X> ddot  (const tensor2d<X> &A, const tensor4 <X> &B);
template<class X> inline          X  ddot  (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline          X  ddot  (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline          X  ddot  (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline          X  ddot  (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline          X  ddot  (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline          X  ddot  (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline          X  ddot  (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline          X  ddot  (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline          X  ddot  (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> dot   (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> dot   (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline vector  <X> dot   (const tensor2 <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> dot   (const tensor2s<X> &A, const vector  <X> &B);
template<class X> inline vector  <X> dot   (const tensor2d<X> &A, const vector  <X> &B);
template<class X> inline vector  <X> dot   (const vector  <X> &A, const tensor2 <X> &B);
template<class X> inline vector  <X> dot   (const vector  <X> &A, const tensor2s<X> &B);
template<class X> inline vector  <X> dot   (const vector  <X> &A, const tensor2d<X> &B);
template<class X> inline          X  dot   (const vector  <X> &A, const vector  <X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> dyadic(const vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> cross (const vector  <X> &A, const vector  <X> &B);

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X> inline tensor2 <X> transpose (const tensor2 <X> &A);
template<class X> inline tensor2s<X> transpose (const tensor2s<X> &A);
template<class X> inline tensor2d<X> transpose (const tensor2d<X> &A);
template<class X> inline tensor4 <X> transpose (const tensor4 <X> &A);
template<class X> inline tensor4 <X> transposeR(const tensor4 <X> &A);
template<class X> inline tensor4 <X> transposeL(const tensor4 <X> &A);
template<class X> inline tensor2 <X> inv       (const tensor2 <X> &A);
template<class X> inline tensor2s<X> inv       (const tensor2s<X> &A);
template<class X> inline tensor2d<X> inv       (const tensor2d<X> &A);
template<class X> inline          X  det       (const tensor2 <X> &A);
template<class X> inline          X  det       (const tensor2s<X> &A);
template<class X> inline          X  det       (const tensor2d<X> &A);
template<class X> inline          X  trace     (const tensor2 <X> &A);
template<class X> inline          X  trace     (const tensor2s<X> &A);
template<class X> inline          X  trace     (const tensor2d<X> &A);

// =================================================================================================
// identity tensors
// =================================================================================================

template<class X> inline tensor4 <X> identity4  (size_t nd);
template<class X> inline tensor4 <X> identity4rt(size_t nd);
template<class X> inline tensor4 <X> identity4II(size_t nd);
template<class X> inline tensor4 <X> identity4s (size_t nd);
template<class X> inline tensor4 <X> identity4d (size_t nd);
template<class X> inline tensor2d<X> identity2  (size_t nd);

// =================================================================================================

}} // namespace ...

#endif

