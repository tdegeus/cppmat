/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR3_H
#define CPPMAT_TENSOR3_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian3d {

// =================================================================================================
// alias name-space with "view" class
// =================================================================================================

namespace map = cppmat::view::cartesian3d;

// =================================================================================================
// cppmat::cartesian3d::tensor4
// =================================================================================================

template<class X>
class tensor4
{
private:

  X m_data[81]; // data container

public:

  // constructor
  tensor4(){};

  // constructor: initialize
  static tensor4<X> Arange();
  static tensor4<X> Zero();
  static tensor4<X> Ones();
  static tensor4<X> Constant(X D);
  static tensor4<X> I();
  static tensor4<X> Irt();
  static tensor4<X> Is();
  static tensor4<X> Id();
  static tensor4<X> II();

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static tensor4<X> Copy(Iterator first, Iterator last);

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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterators
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);
  void setI();
  void setIrt();
  void setIs();
  void setId();
  void setII();

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
// cppmat::cartesian3d::tensor2
// =================================================================================================

template<class X>
class tensor2
{
private:

  X m_data[9];  // data container

public:

  // constructor
  tensor2(){};

  // constructor: initialize
  static tensor2<X> Arange();
  static tensor2<X> Zero();
  static tensor2<X> Ones();
  static tensor2<X> Constant(X D);
  static tensor2<X> I();

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static tensor2<X> Copy(Iterator first, Iterator last);

  // cast into another object
  template<class U> U cast() const;

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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterators
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);
  void setI();

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
  X          inline det   (                    ) const; // determinant
  tensor2<X> inline inv   (                    ) const; // inverse

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

  X m_data[6];  // data container

public:

  // constructor
  tensor2s(){};

  // constructor: initialize
  static tensor2s<X> Arange();
  static tensor2s<X> Zero();
  static tensor2s<X> Ones();
  static tensor2s<X> Constant(X D);
  static tensor2s<X> I();

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static tensor2s<X> Copy(Iterator first, Iterator last);

  // constructor: initialize by copying from external object, stored as "tensor2"
  template<typename Iterator>
  static tensor2s<X> CopyDense(Iterator first, Iterator last);

  // cast into another object
  template<class U> U cast() const;

  // automatic conversion to "tensor2"
  #ifndef CPPMAT_NOCONVERT
  operator tensor2<X> () const;
  #endif

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j);
  const X& operator()(size_t i, size_t j) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterators
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);
  template<typename Iterator> void setCopyDense(Iterator first, Iterator last);
  void setI();

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
  X           inline det   (                    ) const; // determinant
  tensor2s<X> inline inv   (                    ) const; // inverse

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
// cppmat::cartesian3d::tensor2d (diagonal storage of "cppmat::cartesian3d::tensor")
// =================================================================================================

template<class X>
class tensor2d
{
private:

  X m_data[3];  // data container
  X m_zero[1];  // dummy parameter, used to return "0" for any off-diagonal entry

public:

  // constructor
  tensor2d();

  // constructor: initialize
  static tensor2d<X> Arange();
  static tensor2d<X> Zero();
  static tensor2d<X> Ones();
  static tensor2d<X> Constant(X D);
  static tensor2d<X> I();

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static tensor2d<X> Copy(Iterator first, Iterator last);

  // constructor: initialize by copying from external object, stored as "tensor2"
  template<typename Iterator>
  static tensor2d<X> CopyDense(Iterator first, Iterator last);

  // cast into another object
  template<class U> U cast() const;

  // automatic conversion to "tensor2" or "tensor2s"
  #ifndef CPPMAT_NOCONVERT
  operator tensor2 <X> () const;
  operator tensor2s<X> () const;
  #endif

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t i, size_t j);
  const X& operator()(size_t i, size_t j) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterators
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);
  template<typename Iterator> void setCopyDense(Iterator first, Iterator last);
  void setI();

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
  X           inline det   (                    ) const; // determinant
  tensor2d<X> inline inv   (                    ) const; // inverse

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

  X m_data[3];  // data container

public:

  // constructor
  vector(){};

  // constructor: initialize
  static vector<X> Arange();
  static vector<X> Zero();
  static vector<X> Ones();
  static vector<X> Constant(X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static vector<X> Copy(Iterator first, Iterator last);

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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterators
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

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
  vector <X> inline cross (const vector  <X> &B) const; // cross  product

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
// arithmetic operators - mixed class
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
// arithmetic operators - mixed class - view ? normal
// =================================================================================================

template<class X> inline tensor4 <X> operator* (const map::tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator/ (const map::tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator+ (const map::tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator- (const map::tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline tensor4 <X> operator* (const map::tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator/ (const map::tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator+ (const map::tensor4 <X> &A, const          X  &B);
template<class X> inline tensor4 <X> operator- (const map::tensor4 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator* (const map::tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const map::tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator* (const map::tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator/ (const map::tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2 <X> operator* (const map::tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator/ (const map::tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2 <X> &A, const          X  &B);
template<class X> inline tensor2 <X> operator* (const map::tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const map::tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const map::tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const map::tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2s<X> operator* (const map::tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator/ (const map::tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const map::tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2s<X> operator- (const map::tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2s<X> operator* (const map::tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator/ (const map::tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator- (const map::tensor2s<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const map::tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator- (const map::tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2s<X> operator+ (const map::tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const map::tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator+ (const map::tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator- (const map::tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const          X  &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const map::tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline vector  <X> operator* (const map::vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator/ (const map::vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator+ (const map::vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator- (const map::vector  <X> &A, const vector  <X> &B);
template<class X> inline vector  <X> operator* (const map::vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator/ (const map::vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator+ (const map::vector  <X> &A, const          X  &B);
template<class X> inline vector  <X> operator- (const map::vector  <X> &A, const          X  &B);

// =================================================================================================
// arithmetic operators - mixed class - normal ? view
// =================================================================================================

template<class X> inline tensor4 <X> operator* (const tensor4 <X> &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator/ (const tensor4 <X> &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator+ (const tensor4 <X> &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator- (const tensor4 <X> &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator* (const          X  &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator/ (const          X  &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator+ (const          X  &A, const map::tensor4 <X> &B);
template<class X> inline tensor4 <X> operator- (const          X  &A, const map::tensor4 <X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2 <X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2 <X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2 <X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2 <X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2 <X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2 <X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2 <X> operator* (const tensor2s<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const tensor2s<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2s<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2s<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const tensor2d<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const tensor2d<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator* (const          X  &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator/ (const          X  &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator+ (const          X  &A, const map::tensor2 <X> &B);
template<class X> inline tensor2 <X> operator- (const          X  &A, const map::tensor2 <X> &B);
template<class X> inline tensor2s<X> operator* (const tensor2s<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator/ (const tensor2s<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const tensor2s<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2s<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const tensor2s<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2s<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2s<X> operator+ (const tensor2d<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const tensor2d<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator* (const          X  &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator/ (const          X  &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const          X  &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator- (const          X  &A, const map::tensor2s<X> &B);
template<class X> inline tensor2s<X> operator+ (const          X  &A, const map::tensor2d<X> &B);
template<class X> inline tensor2s<X> operator- (const          X  &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator+ (const tensor2d<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator- (const tensor2d<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2d<X> operator/ (const tensor2d<X> &A, const map::tensor2 <X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2d<X> operator/ (const tensor2d<X> &A, const map::tensor2s<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2 <X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const tensor2s<X> &A, const map::tensor2d<X> &B);
template<class X> inline tensor2d<X> operator* (const          X  &A, const map::tensor2d<X> &B);
template<class X> inline vector  <X> operator* (const vector  <X> &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator/ (const vector  <X> &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator+ (const vector  <X> &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator- (const vector  <X> &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator* (const          X  &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator/ (const          X  &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator+ (const          X  &A, const map::vector  <X> &B);
template<class X> inline vector  <X> operator- (const          X  &A, const map::vector  <X> &B);

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

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

