/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TENSOR3_H
#define CPPMAT_VIEW_TENSOR3_H

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
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

template<class X>
class tensor4
{
private:

// pointer to data (points outside)
  const X *m_data;

public:

  // constructors
  tensor4();
  tensor4(const X *D);

  // assignment operator
  tensor4<X>& operator= (const tensor4<X> &D);

  // map external pointer
  void map(const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j, size_t k, size_t l) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // tensor products / operations
  cppmat::cartesian3d::tensor4<X> inline ddot(const tensor4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  cppmat::cartesian3d::tensor2<X> inline ddot(const tensor2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  cppmat::cartesian3d::tensor2<X> inline ddot(const tensor2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  cppmat::cartesian3d::tensor2<X> inline ddot(const tensor2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  cppmat::cartesian3d::tensor4<X> inline T   (                    ) const; // transposition   : C_lkji = A_ijkl
  cppmat::cartesian3d::tensor4<X> inline RT  (                    ) const; // transposition   : C_ijlk = A_ijkl
  cppmat::cartesian3d::tensor4<X> inline LT  (                    ) const; // transposition   : C_jikl = A_ijkl

  // equality operators
  bool operator== (const tensor4<X> &B) const; // A_ijkl == B_ijkl

  // basic algebra
  X norm() const; // sum(| A_ijkl |)

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian3d::tensor2
// =================================================================================================

template<class X>
class tensor2
{
private:

  // pointer to data (points outside)
  const X *m_data;

public:

  // constructors
  tensor2();
  tensor2(const X *D);

  // assignment operator
  tensor2<X>& operator= (const tensor2<X> &D);

  // map external pointer
  void map(const X *D);

  // cast into another object
  template<class U> U cast() const;

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // tensor products / operations
  cppmat::cartesian3d::tensor2<X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2<X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2<X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::vector <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  cppmat::cartesian3d::tensor2<X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X                               inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X                               inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X                               inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  cppmat::cartesian3d::tensor4<X> inline dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4<X> inline dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4<X> inline dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor2<X> inline T     (                    ) const; // transpose       : C_ij   = A_ji
  X                               inline trace (                    ) const; // trace           : A_ii
  X                               inline det   (                    ) const; // determinant
  cppmat::cartesian3d::tensor2<X> inline inv   (                    ) const; // inverse

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
// cppmat::cartesian3d::tensor2s (symmetric storage of "cppmat::cartesian3d::tensor")
// =================================================================================================

template<class X>
class tensor2s
{
private:

  // pointer to data (points outside)
  const X *m_data;

public:

  // constructors
  tensor2s();
  tensor2s(const X *D);

  // assignment operator
  tensor2s<X>& operator= (const tensor2s<X> &D);

  // map external pointer
  void map(const X *D);

  // cast into another object
  template<class U> U cast() const;

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // tensor products / operations
  cppmat::cartesian3d::tensor2 <X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2 <X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2 <X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::vector  <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  cppmat::cartesian3d::tensor2 <X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X                                inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X                                inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X                                inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor2s<X> inline T     (                    ) const; // transpose       : C_ij   = A_ji
  X                                inline trace (                    ) const; // trace           : A_ii
  X                                inline det   (                    ) const; // determinant
  cppmat::cartesian3d::tensor2s<X> inline inv   (                    ) const; // inverse

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
// cppmat::cartesian3d::tensor2d (symmetric storage of "cppmat::cartesian3d::tensor")
// =================================================================================================

template<class X>
class tensor2d
{
private:

  // pointer to data (points outside)
  const X *m_data;
  // dummy parameter, used to return "0" for any off-diagonal entry
  X m_zero[1];

public:

  // constructors
  tensor2d();
  tensor2d(const X *D);

  // assignment operator
  tensor2d<X>& operator= (const tensor2d<X> &D);

  // map external pointer
  void map(const X *D);

  // cast into another object
  template<class U> U cast() const;

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // tensor products / operations
  cppmat::cartesian3d::tensor2 <X> inline dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2 <X> inline dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::tensor2d<X> inline dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  cppmat::cartesian3d::vector  <X> inline dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  cppmat::cartesian3d::tensor2 <X> inline ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X                                inline ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X                                inline ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X                                inline ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2 <X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2s<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor4 <X> inline dyadic(const tensor2d<X> &B) const; // dyadic product : C_ijkl = A_ij * B_kl
  cppmat::cartesian3d::tensor2d<X> inline T     (                    ) const; // transpose      : C_ij   = A_ji
  X                                inline trace (                    ) const; // trace          : A_ii
  X                                inline det   (                    ) const; // determinant
  cppmat::cartesian3d::tensor2d<X> inline inv   (                    ) const; // inverse

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
// cppmat::cartesian3d::vector
// =================================================================================================

template<class X>
class vector
{
private:

  // pointer to data (points outside)
  const X *m_data;

public:

  // constructors
  vector();
  vector(const X *D);

  // assignment operator
  vector<X>& operator= (const vector<X> &D);

  // map external pointer
  void map(const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // tensor products / operations
  X                               inline dot   (const vector  <X> &B) const; // dot    product: C   = A_i*B_i
  cppmat::cartesian3d::vector <X> inline dot   (const tensor2 <X> &B) const; // dot    product: C_j = A_i*B_ij
  cppmat::cartesian3d::vector <X> inline dot   (const tensor2s<X> &B) const; // dot    product: C_j = A_i*B_ij
  cppmat::cartesian3d::vector <X> inline dot   (const tensor2d<X> &B) const; // dot    product: C_j = A_i*B_ij
  cppmat::cartesian3d::tensor2<X> inline dyadic(const vector  <X> &B) const; // dyadic product: C_ij = A_i*B_j
  cppmat::cartesian3d::vector <X> inline cross (const vector  <X> &B) const; // cross  product

  // equality operators
  bool operator== (const vector<X> &B) const; // A_i == B_i

  // basic algebra
  X norm() const;       // sum(| A_i |)
  X length() const;     // sqrt(sum(pow(A_i,2.)))

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X> inline cppmat::cartesian3d::tensor4 <X> operator* (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator/ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator+ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator- (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator* (const tensor4 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator/ (const tensor4 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator+ (const tensor4 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator- (const tensor4 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator* (const          X  &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator/ (const          X  &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator+ (const          X  &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> operator- (const          X  &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator* (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator* (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator* (const tensor2 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator/ (const tensor2 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2 <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator* (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator/ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator* (const          X  &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator/ (const          X  &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator+ (const          X  &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> operator- (const          X  &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator* (const tensor2s<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator/ (const tensor2s<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const tensor2s<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const tensor2s<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const tensor2d<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const tensor2d<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator* (const          X  &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator/ (const          X  &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const          X  &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const          X  &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator+ (const          X  &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2s<X> operator- (const          X  &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2d<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator/ (const tensor2d<X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> operator* (const          X  &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator* (const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator/ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator+ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator- (const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator* (const vector  <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator/ (const vector  <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator+ (const vector  <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator- (const vector  <X> &A, const          X  &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator* (const          X  &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator/ (const          X  &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator+ (const          X  &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> operator- (const          X  &A, const vector  <X> &B);

// =================================================================================================
// tensor products / operations
// =================================================================================================

template<class X> inline cppmat::cartesian3d::tensor4 <X> ddot  (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor2 <X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor2s<X> &A, const tensor4 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> ddot  (const tensor2d<X> &A, const tensor4 <X> &B);
template<class X> inline                               X  ddot  (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline                               X  ddot  (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline                               X  ddot  (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline                               X  ddot  (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline                               X  ddot  (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline                               X  ddot  (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline                               X  ddot  (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline                               X  ddot  (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline                               X  ddot  (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dot   (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2d<X> dot   (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const tensor2 <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const tensor2s<X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const tensor2d<X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const vector  <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const vector  <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> dot   (const vector  <X> &A, const tensor2d<X> &B);
template<class X> inline                               X  dot   (const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline cppmat::cartesian3d::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline cppmat::cartesian3d::tensor2 <X> dyadic(const vector  <X> &A, const vector  <X> &B);
template<class X> inline cppmat::cartesian3d::vector  <X> cross (const vector  <X> &A, const vector  <X> &B);

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X> inline cppmat::cartesian3d::tensor2 <X> transpose (const tensor2 <X> &A);
template<class X> inline cppmat::cartesian3d::tensor2s<X> transpose (const tensor2s<X> &A);
template<class X> inline cppmat::cartesian3d::tensor2d<X> transpose (const tensor2d<X> &A);
template<class X> inline cppmat::cartesian3d::tensor4 <X> transpose (const tensor4 <X> &A);
template<class X> inline cppmat::cartesian3d::tensor4 <X> transposeR(const tensor4 <X> &A);
template<class X> inline cppmat::cartesian3d::tensor4 <X> transposeL(const tensor4 <X> &A);
template<class X> inline cppmat::cartesian3d::tensor2 <X> inv       (const tensor2 <X> &A);
template<class X> inline cppmat::cartesian3d::tensor2s<X> inv       (const tensor2s<X> &A);
template<class X> inline cppmat::cartesian3d::tensor2d<X> inv       (const tensor2d<X> &A);
template<class X> inline                               X  det       (const tensor2 <X> &A);
template<class X> inline                               X  det       (const tensor2s<X> &A);
template<class X> inline                               X  det       (const tensor2d<X> &A);
template<class X> inline                               X  trace     (const tensor2 <X> &A);
template<class X> inline                               X  trace     (const tensor2s<X> &A);
template<class X> inline                               X  trace     (const tensor2d<X> &A);

// =================================================================================================

}}} // namespace ...

#endif

