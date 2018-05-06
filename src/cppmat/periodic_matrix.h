/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX_H
#define CPPMAT_PERIODIC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::periodic::matrix
// =================================================================================================

template<class X>
class matrix
{
private:

  std::vector<X> mData;    // data container
  size_t         mSize=0;  // total size
  size_t         M=0;      // number of rows
  size_t         N=0;      // number of columns
  int            MI=0;     // == M, but int
  int            NI=0;     // == N, but int

public:

  // constructor
  matrix() = default;
  matrix(size_t m, size_t n);

  // constructor: initialize
  static matrix<X> Arange  (size_t m, size_t n);
  static matrix<X> Zero    (size_t m, size_t n);
  static matrix<X> Ones    (size_t m, size_t n);
  static matrix<X> Constant(size_t m, size_t n, X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static matrix<X> Copy(size_t m, size_t n, Iterator first, Iterator last);

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t rows() const;
  size_t cols() const;
  size_t shape(int    i) const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  X&       operator()(int a);
  const X& operator()(int a) const;
  X&       operator()(int a, int b);
  const X& operator()(int a, int b) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of matrix-indices (a,b)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(int a) const;
  size_t compress(int a, int b) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to the first and last entry of a row
  auto beginRow(size_t i);
  auto beginRow(size_t i) const;
  auto endRow(size_t i);
  auto endRow(size_t i) const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using matrix-indices
  auto item(int a);
  auto item(int a) const;
  auto item(int a, int b);
  auto item(int a, int b) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

  // arithmetic operators
  matrix<X>& operator*= (const matrix<X> &B);
  matrix<X>& operator/= (const matrix<X> &B);
  matrix<X>& operator+= (const matrix<X> &B);
  matrix<X>& operator-= (const matrix<X> &B);
  matrix<X>& operator*= (const        X  &B);
  matrix<X>& operator/= (const        X  &B);
  matrix<X>& operator+= (const        X  &B);
  matrix<X>& operator-= (const        X  &B);

  // basic algebra
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  std::vector<size_t> argmin() const;
  std::vector<size_t> argmax() const;
  // - location of the minimum/maximum
  size_t argminIndex() const;
  size_t argmaxIndex() const;
  // - minimum
  X         minCoeff() const;
  vector<X> minCoeff(int    axis) const;
  vector<X> minCoeff(size_t axis) const;
  // - maximum
  X         maxCoeff() const;
  vector<X> maxCoeff(int    axis) const;
  vector<X> maxCoeff(size_t axis) const;
  // - sum
  X         sum() const;
  vector<X> sum(int    axis) const;
  vector<X> sum(size_t axis) const;
  // - mean
  double    mean() const;
  vector<X> mean(int    axis) const;
  vector<X> mean(size_t axis) const;
  // - weighted average
  double    average(const matrix<X> &weights,              bool norm=true) const;
  vector<X> average(const matrix<X> &weights, int    axis, bool norm=true) const;
  vector<X> average(const matrix<X> &weights, size_t axis, bool norm=true) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X> inline matrix<X> operator* (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator- (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator* (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator/ (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator+ (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator- (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator* (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator/ (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator+ (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator- (const        X  &A, const matrix<X> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

