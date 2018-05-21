/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_ARRAY_H
#define CPPMAT_VAR_REGULAR_ARRAY_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// cppmat::array
// =================================================================================================

template<class X>
class array
{
protected:

  static const size_t MAX_DIM=6;    // maximum number of dimensions
  size_t         mSize=0;           // total size == data.size() == prod(shape)
  size_t         mRank=0;           // rank (number of axes)
  size_t         mShape  [MAX_DIM]; // number of entries along each axis
  size_t         mStrides[MAX_DIM]; // stride length for each index
  std::vector<X> mData;             // data container

public:

  // constructor
  array() = default;

  // constructor: allocate, don't initialize
  array(const std::vector<size_t> &shape);

  // constructor: copy
  array(const array<X> &A);

  // named constructor: initialize
  static array<X> Random  (const std::vector<size_t> &shape, X lower=(X)0, X upper=(X)1);
  static array<X> Arange  (const std::vector<size_t> &shape);
  static array<X> Zero    (const std::vector<size_t> &shape);
  static array<X> Ones    (const std::vector<size_t> &shape);
  static array<X> Constant(const std::vector<size_t> &shape, X D);

  // named constructor: copy
  static array<X> Copy(const std::vector<size_t> &shape, const std::vector<X> &D);

  // named constructor: copy
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first);
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first, It last);

  // return plain storage as vector
  std::vector<X> asVector() const;

  // resize
  void resize (const std::vector<size_t> &shape);
  void reshape(const std::vector<size_t> &shape);
  void chrank (size_t rank);

  // get dimensions
  size_t size() const;
  size_t rank() const;
  size_t shape(int    i) const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using array-indices
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
  // N.B. the iterator points to list of array-indices (a,b,c,...)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: array-indices -> plain storage (a,b,c,... -> i)
  size_t compress(size_t a) const;
  size_t compress(size_t a, size_t b) const;
  size_t compress(size_t a, size_t b, size_t c) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d, size_t e) const;
  size_t compress(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const;

  // index operators: plain storage -> array-indices (i -> a,b,c,...)
  std::vector<size_t> decompress(size_t i) const;

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

  // iterator to specific entry: access using array-indices
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

  // initialization
  void setRandom(X lower=(X)0, X upper=(X)1);
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

  // copy to target
  template<typename Iterator> void copyTo(Iterator first) const;
  template<typename Iterator> void copyTo(Iterator first, Iterator last) const;

  // sign change
  array<X> operator- () const;
  array<X> operator+ () const;

  // arithmetic operators
  array<X>& operator*= (const array<X> &B);
  array<X>& operator/= (const array<X> &B);
  array<X>& operator+= (const array<X> &B);
  array<X>& operator-= (const array<X> &B);
  array<X>& operator*= (const       X  &B);
  array<X>& operator/= (const       X  &B);
  array<X>& operator+= (const       X  &B);
  array<X>& operator-= (const       X  &B);

  // absolute value
  array<X> abs() const;

  // norm (sum of absolute values)
  X norm() const;

  // return the indices that would sort an array
  array<size_t> argsort(bool ascending=true) const;

  // location of the minimum/maximum: plain storage (use decompress to convert to indices)
  size_t argmin() const;
  size_t argmax() const;

  // minimum
  X        min() const;
  array<X> min(int    axis) const;
  array<X> min(size_t axis) const;
  array<X> min(const std::vector<int> &axes) const;

  // maximum
  X        max() const;
  array<X> max(int    axis) const;
  array<X> max(size_t axis) const;
  array<X> max(const std::vector<int> &axes) const;

  // sum
  X        sum() const;
  array<X> sum(int    axis) const;
  array<X> sum(size_t axis) const;
  array<X> sum(const std::vector<int> &axes) const;

  // mean
  double   mean() const;
  array<X> mean(int    axis) const;
  array<X> mean(size_t axis) const;
  array<X> mean(const std::vector<int> &axes) const;

  // weighted average
  double   average(const array<X> &weights,                               bool norm=true) const;
  array<X> average(const array<X> &weights, int    axis,                  bool norm=true) const;
  array<X> average(const array<X> &weights, size_t axis,                  bool norm=true) const;
  array<X> average(const array<X> &weights, const std::vector<int> &axes, bool norm=true) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

  // find the plain storage indices of all entries equal to some constant
  std::vector<size_t> where(X D) const;

};

// external arithmetic operators
template<class X> inline array<X> operator* (array<X> A, const array<X> &B);
template<class X> inline array<X> operator/ (array<X> A, const array<X> &B);
template<class X> inline array<X> operator+ (array<X> A, const array<X> &B);
template<class X> inline array<X> operator- (array<X> A, const array<X> &B);
template<class X> inline array<X> operator* (array<X> A, const       X  &B);
template<class X> inline array<X> operator/ (array<X> A, const       X  &B);
template<class X> inline array<X> operator+ (array<X> A, const       X  &B);
template<class X> inline array<X> operator- (array<X> A, const       X  &B);
template<class X> inline array<X> operator* (const X &A,       array<X>  B);
template<class X> inline array<X> operator/ (const X &A, const array<X> &B);
template<class X> inline array<X> operator+ (const X &A,       array<X>  B);
template<class X> inline array<X> operator- (const X &A, const array<X> &B);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

