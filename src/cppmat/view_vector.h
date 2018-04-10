/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_VECTOR_H
#define CPPMAT_VIEW_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

// =================================================================================================
// alias name-space with "normal" class
// =================================================================================================

namespace reg = cppmat::tiny;

// =================================================================================================
// cppmat::view::vector
// =================================================================================================

template<class X, size_t n>
class vector
{
private:

  const X *m_data;  // pointer to data (points outside)
  size_t m_size=n;  // total number of entries

public:

  // constructor
  vector();           // allocate, null-pointer
  vector(const X *D); // allocate, set external pointer

  // map external pointer
  void map(const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(int    i) const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t a) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const vector<X,n> &weights) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t n>
inline reg::vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator* (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline reg::vector<X,n> operator/ (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline reg::vector<X,n> operator+ (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline reg::vector<X,n> operator- (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline reg::vector<X,n> operator* (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator/ (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator+ (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline reg::vector<X,n> operator- (const        X    &A, const vector<X,n> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

