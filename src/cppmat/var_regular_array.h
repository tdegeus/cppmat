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
  size_t         mShape[MAX_DIM];   // number of entries along each axis
  size_t         mStrides[MAX_DIM]; // stride length for each index
  std::vector<X> mData;             // data container
  bool           mPeriodic=false;   // if true: disable bounds-check where possible

public:

  // constructor: default
  array() = default;

  // constructor: allocate, don't initialize
  array(const std::vector<size_t> &shape);

  // constructor: copy from own class (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  array(const cppmat::array<U> &A);

  // constructor: copy from fixed size
  template<size_t rank, size_t i, size_t j, size_t k, size_t l, size_t m, size_t n>
  array(const cppmat::tiny::array<X,rank,i,j,k,l,m,n> &A);

  // constructor: copy from view
  template<size_t rank, size_t i, size_t j, size_t k, size_t l, size_t m, size_t n>
  array(const cppmat::view::array<X,rank,i,j,k,l,m,n> &A);

  // named constructor: initialize
  static array<X> Random  (const std::vector<size_t> &shape, X lower=(X)0, X upper=(X)1);
  static array<X> Arange  (const std::vector<size_t> &shape);
  static array<X> Zero    (const std::vector<size_t> &shape);
  static array<X> Ones    (const std::vector<size_t> &shape);
  static array<X> Constant(const std::vector<size_t> &shape, X D);
  static array<X> Copy    (const std::vector<size_t> &shape, const std::vector<X> &D);

  // named constructor: copy
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first);
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first, It last);

  // return plain storage as vector
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  operator std::vector<U> () const;

  // resize
  void resize (const std::vector<size_t> &shape);
  void resize (const std::vector<size_t> &shape, const X &D);
  void reshape(const std::vector<size_t> &shape);
  void chrank (size_t rank);

  // modify bounds-checks
  void setPeriodic(bool periodic);

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

  // index operators: access using array-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b, T c);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b, T c, T d);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b, T c, T d, T e);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d, T e) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b, T c, T d, T e, T f);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d, T e, T f) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of array-indices (a,b,c,...)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: array-indices -> plain storage (a,b,c,... -> i)
  size_t compress(int a) const;
  size_t compress(int a, int b) const;
  size_t compress(int a, int b, int c) const;
  size_t compress(int a, int b, int c, int d) const;
  size_t compress(int a, int b, int c, int d, int e) const;
  size_t compress(int a, int b, int c, int d, int e, int f) const;

  // index operators: array-indices -> plain storage (a,b,c,... -> i)
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b, T c) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b, T c, T d) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b, T c, T d, T e) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b, T c, T d, T e, T f) const;

  // index operators: plain storage -> array-indices (i -> a,b,c,...)
  std::vector<size_t> decompress(size_t i) const;

  // get index of the midpoint (along a certain axis)
  std::vector<size_t> midpoint() const;
  size_t              midpoint(size_t axis) const;

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

  // iterator to specific entry: access using array-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e, T f);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e, T f) const;

  // slice
  // - allowing negative index
  template<typename T, typename=typename std::enable_if<std::is_integral<T>::value,void>::type>
  array<X> slice(
    const std::vector<T> &a=std::vector<T>(), const std::vector<T> &b=std::vector<T>(),
    const std::vector<T> &c=std::vector<T>(), const std::vector<T> &d=std::vector<T>(),
    const std::vector<T> &e=std::vector<T>(), const std::vector<T> &f=std::vector<T>()
  ) const;
  // - only positive indices
  array<X> slice(
    const std::vector<int> &a=std::vector<int>(), const std::vector<int> &b=std::vector<int>(),
    const std::vector<int> &c=std::vector<int>(), const std::vector<int> &d=std::vector<int>(),
    const std::vector<int> &e=std::vector<int>(), const std::vector<int> &f=std::vector<int>()
  ) const;

  // return padded array
  array<X> pad(const std::vector<size_t> &pad_width, X D=static_cast<X>(0));

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

  // check in an index in within the matrix bound:
  // ( a >= 0 and a < shape[0] ) or ( periodic )
  template<typename T> bool inBounds(T a) const;
  template<typename T> bool inBounds(T a, T b) const;
  template<typename T> bool inBounds(T a, T b, T c) const;
  template<typename T> bool inBounds(T a, T b, T c, T d) const;
  template<typename T> bool inBounds(T a, T b, T c, T d, T e) const;
  template<typename T> bool inBounds(T a, T b, T c, T d, T e, T f) const;

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

  // return array of booleans, based on condition
  array<int> equal        (const       X  &D) const;
  array<int> not_equal    (const       X  &D) const;
  array<int> greater      (const       X  &D) const;
  array<int> greater_equal(const       X  &D) const;
  array<int> less         (const       X  &D) const;
  array<int> less_equal   (const       X  &D) const;
  array<int> equal        (const array<X> &D) const;
  array<int> not_equal    (const array<X> &D) const;
  array<int> greater      (const array<X> &D) const;
  array<int> greater_equal(const array<X> &D) const;
  array<int> less         (const array<X> &D) const;
  array<int> less_equal   (const array<X> &D) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

};

// equality operators
template<class X> bool operator!= (const array<X> &A, const array<X> &B);
template<class X> bool operator== (const array<X> &A, const array<X> &B);

// external arithmetic operators
template<class X> array<X> operator* (const array<X> &A, const array<X> &B);
template<class X> array<X> operator/ (const array<X> &A, const array<X> &B);
template<class X> array<X> operator+ (const array<X> &A, const array<X> &B);
template<class X> array<X> operator- (const array<X> &A, const array<X> &B);
template<class X> array<X> operator* (const array<X> &A, const       X  &B);
template<class X> array<X> operator/ (const array<X> &A, const       X  &B);
template<class X> array<X> operator+ (const array<X> &A, const       X  &B);
template<class X> array<X> operator- (const array<X> &A, const       X  &B);
template<class X> array<X> operator* (const       X  &A, const array<X> &B);
template<class X> array<X> operator/ (const       X  &A, const array<X> &B);
template<class X> array<X> operator+ (const       X  &A, const array<X> &B);
template<class X> array<X> operator- (const       X  &A, const array<X> &B);

// print operator
template<class X> std::ostream& operator<<(std::ostream& out, const array<X>& src);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

