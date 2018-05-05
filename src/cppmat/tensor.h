/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR_H
#define CPPMAT_TENSOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// cppmat::cartesian::tensor4
// =================================================================================================

template<class X>
class tensor4
{
private:

  std::vector<X> mData;    // data container
  size_t         mNd=0;    // number of dimensions (rank == 4 -> shape == [mNd,mNd,mNd,mNd])
  size_t         mSize=0;  // == mData.size() == mNd*mNd*mNd*mNd

public:

  // constructor
  tensor4() = default;
  tensor4(size_t nd);

  // constructor: initialize
  static tensor4<X> Arange  (size_t nd);
  static tensor4<X> Zero    (size_t nd);
  static tensor4<X> Ones    (size_t nd);
  static tensor4<X> Constant(size_t nd, X D);
  static tensor4<X> I       (size_t nd);
  static tensor4<X> Irt     (size_t nd);
  static tensor4<X> Is      (size_t nd);
  static tensor4<X> Id      (size_t nd);
  static tensor4<X> II      (size_t nd);

  // constructor: initialize by copying from external object
  template<typename Iterator> static tensor4<X> Copy(size_t nd, Iterator first, Iterator last);

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

  // index operators: plain storage -> array-indices (i -> a,b,c,d)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b,c,d -> i)
  size_t compress(size_t a, size_t b, size_t c, size_t d) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b, size_t c, size_t d);
  auto item(size_t a, size_t b, size_t c, size_t d) const;

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
  // - tensor norm: sum(| A_ijkl |)
  X norm() const;
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const tensor4<X> &weights, bool norm=true) const;

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

  std::vector<X> mData;    // data container
  size_t         mNd=0;    // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  size_t         mSize=0;  // == mData.size() == mNd*mNd

public:

  // constructor
  tensor2() = default;
  tensor2(size_t nd);

  // constructor: initialize
  static tensor2<X> Arange  (size_t nd);
  static tensor2<X> Zero    (size_t nd);
  static tensor2<X> Ones    (size_t nd);
  static tensor2<X> Constant(size_t nd, X D);
  static tensor2<X> I       (size_t nd);

  // constructor: initialize by copying from external object
  template<typename Iterator> static tensor2<X> Copy(size_t nd, Iterator first, Iterator last);

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

  // index operators: plain storage -> array-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b);
  auto item(size_t a, size_t b) const;

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
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const tensor2<X> &weights, bool norm=true) const;

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

  std::vector<X> mData;    // data container
  size_t         mNd=0;    // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  size_t         mSize=0;  // == mData.size() == (mNd+1)*mNd/2

public:

  // constructor
  tensor2s() = default;
  tensor2s(size_t nd);

  // constructor: initialize
  static tensor2s<X> Arange  (size_t nd);
  static tensor2s<X> Zero    (size_t nd);
  static tensor2s<X> Ones    (size_t nd);
  static tensor2s<X> Constant(size_t nd, X D);
  static tensor2s<X> I       (size_t nd);

  // constructor: initialize by copying from external object
  template<typename Iterator> static tensor2s<X> Copy     (size_t nd, Iterator first, Iterator last);
  template<typename Iterator> static tensor2s<X> CopyDense(size_t nd, Iterator first, Iterator last);

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

  // index operators: plain storage -> array-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b);
  auto item(size_t a, size_t b) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy     (Iterator first, Iterator last);
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
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2s<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // equality operators
  bool operator== (const tensor2 <X> &B) const; // A_ij == A_ij
  bool operator== (const tensor2s<X> &B) const; // A_ij == B_ij
  bool operator== (const tensor2d<X> &B) const; // A_ij == B_ij

  // structure check
  bool isdiagonal() const; // A_ij == 0, i != j

  // basic algebra
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const tensor2s<X> &weights, bool norm=true) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// =================================================================================================
// cppmat::cartesian::tensor2d (diagonal storage of "cppmat::cartesian::tensor")
// =================================================================================================

template<class X>
class tensor2d
{
private:

  std::vector<X> mData;    // data container
  size_t         mNd=0;    // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  size_t         mSize=0;  // == mData.size() == mNd
  X              mZero[1]; // dummy parameter, used to return "0" for any off-diagonal entry

public:

  // constructor
  tensor2d();
  tensor2d(size_t nd);

  // constructor: initialize
  static tensor2d<X> Arange  (size_t nd);
  static tensor2d<X> Zero    (size_t nd);
  static tensor2d<X> Ones    (size_t nd);
  static tensor2d<X> Constant(size_t nd, X D);
  static tensor2d<X> I       (size_t nd);

  // constructor: initialize by copying from external object
  template<typename Iterator> static tensor2d<X> Copy     (size_t nd, Iterator first, Iterator last);
  template<typename Iterator> static tensor2d<X> CopyDense(size_t nd, Iterator first, Iterator last);

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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy     (Iterator first, Iterator last);
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
  X           inline det   (                    ) const; // determinant (only in 2D/3D)
  tensor2d<X> inline inv   (                    ) const; // inverse     (only in 2D/3D)

  // equality operators
  bool operator== (const tensor2 <X> &B) const; // A_ij == B_ij
  bool operator== (const tensor2s<X> &B) const; // A_ij == A_ij
  bool operator== (const tensor2d<X> &B) const; // A_ij == B_ij

  // basic algebra
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const tensor2d<X> &weights, bool norm=true) const;

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

  std::vector<X> mData;    // data container
  size_t         mNd=0;    // number of dimensions (rank == 1 -> shape == [mNd])
  size_t         mSize=0;  // == mData.size() == mNd

public:

  // constructor
  vector() = default;
  vector(size_t nd);

  // constructor: initialize
  static vector<X> Arange  (size_t nd);
  static vector<X> Zero    (size_t nd);
  static vector<X> Ones    (size_t nd);
  static vector<X> Constant(size_t nd, X D);

  // constructor: initialize by copying from external object
  template<typename Iterator> static vector<X> Copy(size_t nd, Iterator first, Iterator last);

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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a);
  auto item(size_t a) const;

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
  vector <X> inline cross (const vector  <X> &B) const; // cross  product (only in 3D)

  // equality operators
  bool operator== (const vector<X> &B) const; // A_i == B_i

  // basic algebra
  // - vector norm: sum(| A_ij |)
  X norm() const;
  // - vector length: sqrt(sum(pow(A_i,2.)))
  X length() const;
  // - set vector length to unity: A_i /= A.length()
  void setUnitLength();
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const vector<X> &weights, bool norm=true) const;

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

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

