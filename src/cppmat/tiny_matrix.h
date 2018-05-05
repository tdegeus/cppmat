/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_MATRIX2_H
#define CPPMAT_TINY_MATRIX2_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix
{
private:

  X mData[M*N];                  // data container
  static const size_t mSize=M*N; // total size

public:

  // constructor
  matrix() = default;

  // constructor: initialize
  static matrix<X,M,N> Arange();
  static matrix<X,M,N> Zero();
  static matrix<X,M,N> Ones();
  static matrix<X,M,N> Constant(X D);

  // constructor: initialize by copying from external object
  template<typename Iterator> static matrix<X,M,N> Copy(Iterator first);
  template<typename Iterator> static matrix<X,M,N> Copy(Iterator first, Iterator last);

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
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of matrix-indices (a,b)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(size_t a) const;
  size_t compress(size_t a, size_t b) const;

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
  auto item(size_t a);
  auto item(size_t a) const;
  auto item(size_t a, size_t b);
  auto item(size_t a, size_t b) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

  // arithmetic operators
  matrix<X,M,N>& operator*= (const matrix<X,M,N> &B);
  matrix<X,M,N>& operator/= (const matrix<X,M,N> &B);
  matrix<X,M,N>& operator+= (const matrix<X,M,N> &B);
  matrix<X,M,N>& operator-= (const matrix<X,M,N> &B);
  matrix<X,M,N>& operator*= (const        X      &B);
  matrix<X,M,N>& operator/= (const        X      &B);
  matrix<X,M,N>& operator+= (const        X      &B);
  matrix<X,M,N>& operator-= (const        X      &B);

  // basic algebra
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
  double average(const matrix<X,M,N> &weights, bool norm=true) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator* (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator- (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator* (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator- (const matrix<X,M,N> &A, const        X      &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator* (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator/ (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator+ (const        X      &A, const matrix<X,M,N> &B);

template<class X, size_t M, size_t N>
inline matrix<X,M,N> operator- (const        X      &A, const matrix<X,M,N> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

