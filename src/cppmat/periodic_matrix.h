/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX_H
#define CPPMAT_PERIODIC_MATRIX_H

#include "macros.h"

namespace cppmat {
namespace periodic {

#define MAX_DIM 6

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template<class X>
class matrix
{
private:

  std::vector<X> m_data;             // data container
  size_t         m_ndim=0;           // actual number of dimensions
  size_t         m_size=0;           // total number of entries == data.size() == prod(shape)
  size_t         m_shape[MAX_DIM];   // number of entries in each dimensions
  int            m_shape_i[MAX_DIM]; // == m_shape, but with int's
  size_t         m_strides[MAX_DIM]; // stride length for each index

public:

  // constructors
  matrix(){};
  matrix(const std::vector<size_t> &shape);
  matrix(const std::vector<size_t> &shape, X D);

  // resize
  void resize (const std::vector<size_t> &shape);
  void reshape(const std::vector<size_t> &shape);
  void chdim  (size_t ndim);

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
  X&       operator()(int a);
  const X& operator()(int a) const;
  X&       operator()(int a, int b);
  const X& operator()(int a, int b) const;
  X&       operator()(int a, int b, int c);
  const X& operator()(int a, int b, int c) const;
  X&       operator()(int a, int b, int c, int d);
  const X& operator()(int a, int b, int c, int d) const;
  X&       operator()(int a, int b, int c, int d, int e);
  const X& operator()(int a, int b, int c, int d, int e) const;
  X&       operator()(int a, int b, int c, int d, int e, int f);
  const X& operator()(int a, int b, int c, int d, int e, int f) const;

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
  matrix<X>& operator*= (const matrix<X> &B);
  matrix<X>& operator/= (const matrix<X> &B);
  matrix<X>& operator+= (const matrix<X> &B);
  matrix<X>& operator-= (const matrix<X> &B);
  matrix<X>& operator*= (const        X  &B);
  matrix<X>& operator/= (const        X  &B);
  matrix<X>& operator+= (const        X  &B);
  matrix<X>& operator-= (const        X  &B);

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const matrix<X> &weights) const;

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

#endif

