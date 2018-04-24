/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX2_H
#define CPPMAT_PERIODIC_MATRIX2_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::periodic::matrix2
// =================================================================================================

template<class X>
class matrix2
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_m=0;     // number of rows
  size_t         m_n=0;     // number of columns
  int            m_m_i=0;   // == m_m, but int
  int            m_n_i=0;   // == m_n, but int
  size_t         m_size=0;  // total size

public:

  // constructor
  matrix2(){};
  matrix2(size_t m, size_t n);

  // constructor: initialize
  static matrix2<X> Arange  (size_t m, size_t n);
  static matrix2<X> Zero    (size_t m, size_t n);
  static matrix2<X> Ones    (size_t m, size_t n);
  static matrix2<X> Constant(size_t m, size_t n, X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static matrix2<X> Copy(size_t m, size_t n, Iterator first, Iterator last);

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

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
  X&       operator()(int a);
  const X& operator()(int a) const;
  X&       operator()(int a, int b);
  const X& operator()(int a, int b) const;

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
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

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
  X      minCoeff() const;
  X      maxCoeff() const;
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

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

