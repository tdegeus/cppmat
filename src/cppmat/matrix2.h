/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MATRIX2_H
#define CPPMAT_MATRIX2_H

#include "macros.h"

namespace cppmat {

// =================================================================================================
// cppmat::matrix2
// =================================================================================================

template<class X>
class matrix2
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_m=0;     // number of rows
  size_t         m_n=0;     // number of columns
  size_t         m_size=0;  // total size

public:

  // constructors
  matrix2(){};
  matrix2(size_t m, size_t n);
  matrix2(size_t m, size_t n, X D);

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

  // arithmetic operators
  matrix2<X>& operator*= (const matrix2<X> &B);
  matrix2<X>& operator/= (const matrix2<X> &B);
  matrix2<X>& operator+= (const matrix2<X> &B);
  matrix2<X>& operator-= (const matrix2<X> &B);
  matrix2<X>& operator*= (const         X  &B);
  matrix2<X>& operator/= (const         X  &B);
  matrix2<X>& operator+= (const         X  &B);
  matrix2<X>& operator-= (const         X  &B);

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const matrix2<X> &weights) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X> inline matrix2<X> operator* (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator/ (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator+ (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator- (const matrix2<X> &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator* (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator/ (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator+ (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator- (const matrix2<X> &A, const         X  &B);
template<class X> inline matrix2<X> operator* (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator/ (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator+ (const         X  &A, const matrix2<X> &B);
template<class X> inline matrix2<X> operator- (const         X  &A, const matrix2<X> &B);

// =================================================================================================

} // namespace ...

#endif

