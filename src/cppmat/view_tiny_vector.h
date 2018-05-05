/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TINY_VECTOR_H
#define CPPMAT_VIEW_TINY_VECTOR_H

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
// cppmat::view::vector
// =================================================================================================

template<class X, size_t N>
class vector
{
private:

  const X *mData=nullptr;      // pointer to data (points outside)
  static const size_t mSize=N; // total size

public:

  // constructor
  vector() = default;

  // constructor: map external pointer
  static vector<X,N> Map(const X *D);

  // reset external pointer
  void setMap(const X *D);

  // information without constructing
  static size_t Size();

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(int    i) const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using list-index
  const X& operator()(size_t a) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using list-index
  auto item(size_t a) const;
  // basic algebra
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  size_t argmin() const;
  size_t argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const vector<X,N> &weights, bool norm=true) const;

  // find all non-zero entries
  cppmat::vector<size_t> where() const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t N>
inline reg::vector<X,N> operator* (const vector<X,N> &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const vector<X,N> &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const vector<X,N> &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator- (const vector<X,N> &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator* (const vector<X,N> &A, const        X    &B);

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const vector<X,N> &A, const        X    &B);

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const vector<X,N> &A, const        X    &B);

template<class X, size_t N>
inline reg::vector<X,N> operator- (const vector<X,N> &A, const        X    &B);

template<class X, size_t N>
inline reg::vector<X,N> operator* (const        X    &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const        X    &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const        X    &A, const vector<X,N> &B);

template<class X, size_t N>
inline reg::vector<X,N> operator- (const        X    &A, const vector<X,N> &B);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

