/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_MATRIX2_H
#define CPPMAT_VIEW_MATRIX2_H

#include "cppmat.h"

namespace cppmat {
namespace view {

// =================================================================================================
// cppmat::view::matrix2
// =================================================================================================

template<class X, size_t m, size_t n>
class matrix2
{
private:

  // pointer to data (points outside)
  const X *m_data;
  // total number of entries
  size_t m_size=m*n;

public:

  // constructors
  matrix2();
  matrix2(const X *D);

  // assignment operator
  matrix2<X,m,n>& operator= (const matrix2<X,m,n> &D);

  // map external pointer
  void map(const X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  const X& operator()(size_t a) const;
  const X& operator()(size_t a, size_t b) const;

  // pointer / iterators
  const X* data() const;
  auto     begin() const;
  auto     end() const;

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const matrix2<X,m,n> &weights) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
// TODO for tiny as well
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const         X      &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator* (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator/ (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator+ (const         X      &A, const matrix2<X,m,n> &B);
template<class X, size_t m, size_t n> inline cppmat::tiny::matrix2<X,m,n> operator- (const         X      &A, const matrix2<X,m,n> &B);

// =================================================================================================

}} // namespace ...

#endif

