/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TINY_MATRIX_H
#define CPPMAT_VIEW_TINY_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace tiny {

// =================================================================================================
// alias name-space with "normal" class
// =================================================================================================

namespace reg = cppmat::tiny;

// =================================================================================================
// cppmat::view::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix
{
private:

  const X *mData=nullptr;        // pointer to data (points outside)
  static const size_t mSize=M*N; // total size

public:

  // constructor
  matrix() = default;

  // constructor: map external pointer
  static matrix<X,M,N> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

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
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  const X& operator()(size_t a) const;
  const X& operator()(size_t a, size_t b) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of matrix-indices (a,b)
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(size_t a) const;
  size_t compress(size_t a, size_t b) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to the first and last entry of a row
  auto beginRow(size_t i) const;
  auto endRow(size_t i) const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using matrix-indices
  auto item(size_t a) const;
  auto item(size_t a, size_t b) const;

  // basic algebra
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
  double average(const matrix<X,M,N> &weights, bool norm=true) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator* (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator/ (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator+ (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline reg::matrix<X,M,N> operator- (const        X      &A, const matrix<X,M,N> &B);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

