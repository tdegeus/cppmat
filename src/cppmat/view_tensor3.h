/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TENSOR3_H
#define CPPMAT_VIEW_TENSOR3_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian3d {

// =================================================================================================
// alias name-space with "normal" class
// =================================================================================================

namespace reg = cppmat::cartesian3d;

// =================================================================================================
// cppmat::view::cartesian3d::tensor4
// =================================================================================================

template<class X>
class tensor4
{
private:

  const X *mData=nullptr;       // pointer to data (points outside)
  static const size_t mNd=3;    // number of dimensions (rank == 4 -> shape == [mNd,mNd,mNd,mNd])
  static const size_t mSize=81; // == mData.size() == mNd*mNd*mNd*mNd

public:

  // constructor
  tensor4() = default;

  // constructor: map external pointer
  static tensor4<X> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j, size_t k, size_t l) const;

  // index operators: plain storage -> array-indices (i -> a,b,c,d)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b,c,d -> i)
  size_t compress(size_t a, size_t b, size_t c, size_t d) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b, size_t c, size_t d) const;
  // tensor products / operations
  reg::tensor4<X> inline ddot(const tensor4 <X> &B) const;
  reg::tensor2<X> inline ddot(const tensor2 <X> &B) const;
  reg::tensor2<X> inline ddot(const tensor2s<X> &B) const;
  reg::tensor2<X> inline ddot(const tensor2d<X> &B) const;
  reg::tensor4<X> inline T   (                    ) const;
  reg::tensor4<X> inline RT  (                    ) const;
  reg::tensor4<X> inline LT  (                    ) const;

  // equality operators
  bool operator== (const      tensor4<X> &B) const;
  bool operator== (const reg::tensor4<X> &B) const;

  // basic algebra
  // - tensor norm: sum(| A_ijkl |)
  X norm() const;
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
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
// cppmat::view::cartesian3d::tensor2
// =================================================================================================

template<class X>
class tensor2
{
private:

  const X *mData=nullptr;      // pointer to data (points outside)
  static const size_t mNd=3;   // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  static const size_t mSize=9; // == mData.size() == mNd*mNd

public:

  // constructor
  tensor2() = default;

  // constructor: map external pointer
  static tensor2<X> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

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

  // index operators: plain storage -> array-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b) const;

  // tensor products / operations
  reg::tensor2<X> inline dot   (const tensor2 <X> &B) const;
  reg::tensor2<X> inline dot   (const tensor2s<X> &B) const;
  reg::tensor2<X> inline dot   (const tensor2d<X> &B) const;
  reg::vector <X> inline dot   (const vector  <X> &B) const;
  reg::tensor2<X> inline ddot  (const tensor4 <X> &B) const;
  X               inline ddot  (const tensor2 <X> &B) const;
  X               inline ddot  (const tensor2s<X> &B) const;
  X               inline ddot  (const tensor2d<X> &B) const;
  reg::tensor4<X> inline dyadic(const tensor2 <X> &B) const;
  reg::tensor4<X> inline dyadic(const tensor2s<X> &B) const;
  reg::tensor4<X> inline dyadic(const tensor2d<X> &B) const;
  reg::tensor2<X> inline T     (                    ) const;
  X               inline trace (                    ) const;
  X               inline det   (                    ) const;
  reg::tensor2<X> inline inv   (                    ) const;

  // equality operators
  bool operator== (const      tensor2 <X> &B) const;
  bool operator== (const      tensor2s<X> &B) const;
  bool operator== (const      tensor2d<X> &B) const;
  bool operator== (const reg::tensor2 <X> &B) const;
  bool operator== (const reg::tensor2s<X> &B) const;
  bool operator== (const reg::tensor2d<X> &B) const;

  // structure check
  bool issymmetric() const;
  bool isdiagonal() const;

  // basic algebra
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
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
// cppmat::view::cartesian3d::tensor2s
// =================================================================================================

template<class X>
class tensor2s
{
private:

  const X *mData=nullptr;      // pointer to data (points outside)
  static const size_t mNd=3;   // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  static const size_t mSize=6; // == mData.size() == (mNd+1)*mNd/2

public:

  // constructor
  tensor2s() = default;

  // constructor: map external pointer
  static tensor2s<X> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

  // cast into another object
  template<class U> U cast() const;

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j) const;

  // index operators: plain storage -> array-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: array-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a, size_t b) const;

  // tensor products / operations
  reg::tensor2 <X> inline dot   (const tensor2 <X> &B) const;
  reg::tensor2 <X> inline dot   (const tensor2s<X> &B) const;
  reg::tensor2 <X> inline dot   (const tensor2d<X> &B) const;
  reg::vector  <X> inline dot   (const vector  <X> &B) const;
  reg::tensor2 <X> inline ddot  (const tensor4 <X> &B) const;
  X                inline ddot  (const tensor2 <X> &B) const;
  X                inline ddot  (const tensor2s<X> &B) const;
  X                inline ddot  (const tensor2d<X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2 <X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2s<X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2d<X> &B) const;
  reg::tensor2s<X> inline T     (                    ) const;
  X                inline trace (                    ) const;
  X                inline det   (                    ) const;
  reg::tensor2s<X> inline inv   (                    ) const;

  // equality operators
  bool operator== (const      tensor2 <X> &B) const;
  bool operator== (const      tensor2s<X> &B) const;
  bool operator== (const      tensor2d<X> &B) const;
  bool operator== (const reg::tensor2 <X> &B) const;
  bool operator== (const reg::tensor2s<X> &B) const;
  bool operator== (const reg::tensor2d<X> &B) const;

  // structure check
  bool isdiagonal() const; // A_ij == 0, i != j

  // basic algebra
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
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
// cppmat::view::cartesian3d::tensor2d
// =================================================================================================

template<class X>
class tensor2d
{
private:

  const X *mData=nullptr;      // pointer to data (points outside)
  static const size_t mNd=3;   // number of dimensions (rank == 2 -> shape == [mNd,mNd])
  static const size_t mSize=3; // == mData.size() == mNd
  X mZero[1];                  // dummy parameter, used to return "0" for any off-diagonal entry

public:

  // constructor
  tensor2d();

  // constructor: map external pointer
  static tensor2d<X> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

  // cast into another object
  template<class U> U cast() const;

  // get dimensions
  size_t size() const;
  size_t ndim() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i, size_t j) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // tensor products / operations
  reg::tensor2 <X> inline dot   (const tensor2 <X> &B) const;
  reg::tensor2 <X> inline dot   (const tensor2s<X> &B) const;
  reg::tensor2d<X> inline dot   (const tensor2d<X> &B) const;
  reg::vector  <X> inline dot   (const vector  <X> &B) const;
  reg::tensor2 <X> inline ddot  (const tensor4 <X> &B) const;
  X                inline ddot  (const tensor2 <X> &B) const;
  X                inline ddot  (const tensor2s<X> &B) const;
  X                inline ddot  (const tensor2d<X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2 <X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2s<X> &B) const;
  reg::tensor4 <X> inline dyadic(const tensor2d<X> &B) const;
  reg::tensor2d<X> inline T     (                    ) const;
  X                inline trace (                    ) const;
  X                inline det   (                    ) const;
  reg::tensor2d<X> inline inv   (                    ) const;

  // equality operators
  bool operator== (const      tensor2 <X> &B) const;
  bool operator== (const      tensor2s<X> &B) const;
  bool operator== (const      tensor2d<X> &B) const;
  bool operator== (const reg::tensor2 <X> &B) const;
  bool operator== (const reg::tensor2s<X> &B) const;
  bool operator== (const reg::tensor2d<X> &B) const;

  // basic algebra
  // - tensor norm: sum(| A_ij |)
  X norm() const;
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
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
// cppmat::view::cartesian3d::vector
// =================================================================================================

template<class X>
class vector
{
private:

  const X *mData=nullptr;      // pointer to data (points outside)
  static const size_t mNd=3;   // number of dimensions (rank == 1 -> shape == [mNd])
  static const size_t mSize=3; // == mData.size() == mNd

public:

  // constructor
  vector() = default;

  // constructor: map external pointer
  static vector<X> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t i) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(size_t a) const;

  // tensor products / operations
  X               inline dot   (const vector  <X> &B) const;
  reg::vector <X> inline dot   (const tensor2 <X> &B) const;
  reg::vector <X> inline dot   (const tensor2s<X> &B) const;
  reg::vector <X> inline dot   (const tensor2d<X> &B) const;
  reg::tensor2<X> inline dyadic(const vector  <X> &B) const;
  reg::vector <X> inline cross (const vector  <X> &B) const;

  // equality operators
  bool operator== (const      vector<X> &B) const;
  bool operator== (const reg::vector<X> &B) const;

  // basic algebra
  // - vector norm: sum(| A_ij |)
  X norm() const;
  // - vector length: sqrt(sum(pow(A_i,2.)))
  X length() const;
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
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

template<class X> inline reg::tensor4 <X> operator* (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator/ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator+ (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator- (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator* (const tensor4 <X> &A, const          X  &B);
template<class X> inline reg::tensor4 <X> operator/ (const tensor4 <X> &A, const          X  &B);
template<class X> inline reg::tensor4 <X> operator+ (const tensor4 <X> &A, const          X  &B);
template<class X> inline reg::tensor4 <X> operator- (const tensor4 <X> &A, const          X  &B);
template<class X> inline reg::tensor4 <X> operator* (const          X  &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator/ (const          X  &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator+ (const          X  &A, const tensor4 <X> &B);
template<class X> inline reg::tensor4 <X> operator- (const          X  &A, const tensor4 <X> &B);
template<class X> inline reg::tensor2 <X> operator* (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator* (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> operator/ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> operator* (const tensor2 <X> &A, const          X  &B);
template<class X> inline reg::tensor2 <X> operator/ (const tensor2 <X> &A, const          X  &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2 <X> &A, const          X  &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2 <X> &A, const          X  &B);
template<class X> inline reg::tensor2 <X> operator* (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator/ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator- (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator* (const          X  &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator/ (const          X  &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator+ (const          X  &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> operator- (const          X  &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2s<X> operator* (const tensor2s<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator/ (const tensor2s<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator+ (const tensor2d<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator- (const tensor2d<X> &A, const          X  &B);
template<class X> inline reg::tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator* (const          X  &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator/ (const          X  &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator+ (const          X  &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator- (const          X  &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2s<X> operator+ (const          X  &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2s<X> operator- (const          X  &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const          X  &B);
template<class X> inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const          X  &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2d<X> operator* (const          X  &A, const tensor2d<X> &B);
template<class X> inline reg::vector  <X> operator* (const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator/ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator+ (const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator- (const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator* (const vector  <X> &A, const          X  &B);
template<class X> inline reg::vector  <X> operator/ (const vector  <X> &A, const          X  &B);
template<class X> inline reg::vector  <X> operator+ (const vector  <X> &A, const          X  &B);
template<class X> inline reg::vector  <X> operator- (const vector  <X> &A, const          X  &B);
template<class X> inline reg::vector  <X> operator* (const          X  &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator/ (const          X  &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator+ (const          X  &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> operator- (const          X  &A, const vector  <X> &B);

// =================================================================================================
// tensor products / operations
// =================================================================================================

template<class X> inline reg::tensor4 <X> ddot  (const tensor4 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor4 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor2 <X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor2s<X> &A, const tensor4 <X> &B);
template<class X> inline reg::tensor2 <X> ddot  (const tensor2d<X> &A, const tensor4 <X> &B);
template<class X> inline               X  ddot  (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline               X  ddot  (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline               X  ddot  (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline               X  ddot  (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline               X  ddot  (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline               X  ddot  (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline               X  ddot  (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline               X  ddot  (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline               X  ddot  (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor2 <X> dot   (const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor2d<X> dot   (const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::vector  <X> dot   (const tensor2 <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> dot   (const tensor2s<X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> dot   (const tensor2d<X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> dot   (const vector  <X> &A, const tensor2 <X> &B);
template<class X> inline reg::vector  <X> dot   (const vector  <X> &A, const tensor2s<X> &B);
template<class X> inline reg::vector  <X> dot   (const vector  <X> &A, const tensor2d<X> &B);
template<class X> inline               X  dot   (const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2 <X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2 <X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B);
template<class X> inline reg::tensor4 <X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B);
template<class X> inline reg::tensor2 <X> dyadic(const vector  <X> &A, const vector  <X> &B);
template<class X> inline reg::vector  <X> cross (const vector  <X> &A, const vector  <X> &B);

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X> inline reg::tensor2 <X> transpose (const tensor2 <X> &A);
template<class X> inline reg::tensor2s<X> transpose (const tensor2s<X> &A);
template<class X> inline reg::tensor2d<X> transpose (const tensor2d<X> &A);
template<class X> inline reg::tensor4 <X> transpose (const tensor4 <X> &A);
template<class X> inline reg::tensor4 <X> transposeR(const tensor4 <X> &A);
template<class X> inline reg::tensor4 <X> transposeL(const tensor4 <X> &A);
template<class X> inline reg::tensor2 <X> inv       (const tensor2 <X> &A);
template<class X> inline reg::tensor2s<X> inv       (const tensor2s<X> &A);
template<class X> inline reg::tensor2d<X> inv       (const tensor2d<X> &A);
template<class X> inline               X  det       (const tensor2 <X> &A);
template<class X> inline               X  det       (const tensor2s<X> &A);
template<class X> inline               X  det       (const tensor2d<X> &A);
template<class X> inline               X  trace     (const tensor2 <X> &A);
template<class X> inline               X  trace     (const tensor2s<X> &A);
template<class X> inline               X  trace     (const tensor2d<X> &A);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

