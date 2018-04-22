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
// cppmat::tiny::matrix2
// =================================================================================================

template<class X, size_t m, size_t n>
class matrix2
{
private:

  X      m_data[m*n]; // data container
  size_t m_size=m*n;  // total number of entries

public:

  // constructor
  matrix2(){};

  // constructor: initialize
  static matrix2<X,m,n> Arange();
  static matrix2<X,m,n> Zero();
  static matrix2<X,m,n> Ones();
  static matrix2<X,m,n> Constant(X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static matrix2<X,m,n> Copy(Iterator first, Iterator last);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(int    i) const;
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

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);

  // arithmetic operators
  matrix2<X,m,n>& operator*= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator/= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator+= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator-= (const matrix2<X,m,n> &B);
  matrix2<X,m,n>& operator*= (const         X      &B);
  matrix2<X,m,n>& operator/= (const         X      &B);
  matrix2<X,m,n>& operator+= (const         X      &B);
  matrix2<X,m,n>& operator-= (const         X      &B);

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
template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const matrix2<X,m,n> &A, const         X      &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const matrix2<X,m,n> &A, const         X      &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const matrix2<X,m,n> &A, const         X      &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const matrix2<X,m,n> &A, const         X      &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator* (const         X      &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator/ (const         X      &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator+ (const         X      &A, const matrix2<X,m,n> &B);

template<class X, size_t m, size_t n>
inline matrix2<X,m,n> operator- (const         X      &A, const matrix2<X,m,n> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

