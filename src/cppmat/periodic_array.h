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
// cppmat::matrix
// =================================================================================================

template<class X>
class matrix
{
private:

  static const size_t MAX_DIM=6;     // maximum number of dimensions
  std::vector<X> m_data;             // data container
  size_t         m_ndim=0;           // actual number of dimensions
  size_t         m_size=0;           // total number of entries == data.size() == prod(shape)
  size_t         m_shape[MAX_DIM];   // number of entries in each dimensions
  int            m_shape_i[MAX_DIM]; // == m_shape, but with int's
  size_t         m_strides[MAX_DIM]; // stride length for each index

public:

  // constructor
  matrix(){};
  matrix(const std::vector<size_t> &shape);

  // constructor: initialize
  static matrix<X> Arange  (const std::vector<size_t> &shape);
  static matrix<X> Zero    (const std::vector<size_t> &shape);
  static matrix<X> Ones    (const std::vector<size_t> &shape);
  static matrix<X> Constant(const std::vector<size_t> &shape, X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static matrix<X> Copy(const std::vector<size_t> &shape, Iterator first, Iterator last);

  // resize
  void resize (const std::vector<size_t> &shape);
  void reshape(const std::vector<size_t> &shape);
  void chdim  (size_t ndim);

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
  X&       operator()(int a, int b, int c);
  const X& operator()(int a, int b, int c) const;
  X&       operator()(int a, int b, int c, int d);
  const X& operator()(int a, int b, int c, int d) const;
  X&       operator()(int a, int b, int c, int d, int e);
  const X& operator()(int a, int b, int c, int d, int e) const;
  X&       operator()(int a, int b, int c, int d, int e, int f);
  const X& operator()(int a, int b, int c, int d, int e, int f) const;

  // index operators: access using iterator
  // N.B. the iterator points to a list of indices (a,b,c,...)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: plain storage -> matrix indices (i -> a,b,c,...)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: matrix indices -> plain storage (a,b,c,... -> i)
  size_t compress(int a) const;
  size_t compress(int a, int b) const;
  size_t compress(int a, int b, int c) const;
  size_t compress(int a, int b, int c, int d) const;
  size_t compress(int a, int b, int c, int d, int e) const;
  size_t compress(int a, int b, int c, int d, int e, int f) const;

  // pointer to data
  X*       data();
  const X* data() const;

  // iterator to first and last entry
  auto begin();
  auto begin() const;
  auto end();
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i);
  auto index(size_t i) const;

  // iterator to specific entry: access using matrix indices
  auto item(int a);
  auto item(int a) const;
  auto item(int a, int b);
  auto item(int a, int b) const;
  auto item(int a, int b, int c);
  auto item(int a, int b, int c) const;
  auto item(int a, int b, int c, int d);
  auto item(int a, int b, int c, int d) const;
  auto item(int a, int b, int c, int d, int e);
  auto item(int a, int b, int c, int d, int e) const;
  auto item(int a, int b, int c, int d, int e, int f);
  auto item(int a, int b, int c, int d, int e, int f) const;

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
  // - minimum/maximum
  X         minCoeff() const;
  X         maxCoeff() const;
  // - sum
  X         sum() const;
  matrix<X> sum(int    axis) const;
  matrix<X> sum(size_t axis) const;
  matrix<X> sum(const std::vector<int> &axes) const;
  // - mean
  double    mean() const;
  matrix<X> mean(int    axis) const;
  matrix<X> mean(size_t axis) const;
  matrix<X> mean(const std::vector<int> &axes) const;
  // - weighted average
  double    average(const matrix<X> &weights,                               bool norm=true) const;
  matrix<X> average(const matrix<X> &weights, int    axis,                  bool norm=true) const;
  matrix<X> average(const matrix<X> &weights, size_t axis,                  bool norm=true) const;
  matrix<X> average(const matrix<X> &weights, const std::vector<int> &axes, bool norm=true) const;

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

