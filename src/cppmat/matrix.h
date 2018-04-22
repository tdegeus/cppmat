/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MATRIX_H
#define CPPMAT_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

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
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;
  X&       operator()(size_t a, size_t b, size_t c);
  const X& operator()(size_t a, size_t b, size_t c) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d);
  const X& operator()(size_t a, size_t b, size_t c, size_t d) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d, size_t e);
  const X& operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const;
  X&       operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
  const X& operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const;

  // index operators: access using iterator
  // N.B. the iterator points to a list of indices (a,b,...)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: plain storage -> matrix indices (i -> a,b,c,...)
  std::vector<size_t> decompress(size_t i) const;

  // index operators: matrix indices -> plain storage (a,b,c,... -> i)
  size_t compress(size_t a) const;
  size_t compress(size_t a, size_t b) const;
  size_t compress(size_t a, size_t b, size_t c) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d, size_t e) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const;

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
  auto item(size_t a);
  auto item(size_t a) const;
  auto item(size_t a, size_t b);
  auto item(size_t a, size_t b) const;
  auto item(size_t a, size_t b, size_t c);
  auto item(size_t a, size_t b, size_t c) const;
  auto item(size_t a, size_t b, size_t c, size_t d);
  auto item(size_t a, size_t b, size_t c, size_t d) const;
  auto item(size_t a, size_t b, size_t c, size_t d, size_t e);
  auto item(size_t a, size_t b, size_t c, size_t d, size_t e) const;
  auto item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f);
  auto item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);

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
  X         min() const;
  X         max() const;
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
  double    average(const matrix<X> &weights) const;
  matrix<X> average(const matrix<X> &weights, int    axis) const;
  matrix<X> average(const matrix<X> &weights, size_t axis) const;
  matrix<X> average(const matrix<X> &weights, const std::vector<int> &axes) const;

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

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

