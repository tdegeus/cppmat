/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_ARRAY_H
#define CPPMAT_FIX_REGULAR_ARRAY_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::array
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
class array
{
  static_assert( RANK >= 1 or I == 1, "Insufficient rank" );
  static_assert( RANK >= 2 or J == 1, "Insufficient rank" );
  static_assert( RANK >= 3 or K == 1, "Insufficient rank" );
  static_assert( RANK >= 4 or L == 1, "Insufficient rank" );
  static_assert( RANK >= 5 or M == 1, "Insufficient rank" );
  static_assert( RANK >= 6 or N == 1, "Insufficient rank" );

protected:

  static const size_t MAX_DIM=6;          // maximum number of dimensions
  static const size_t mSize=I*J*K*L*M*N;  // total size == data.size() == prod(shape)
  static const size_t mRank=RANK;         // rank (number of axes)
  size_t              mShape  [MAX_DIM];  // number of entries along each axis
  size_t              mStrides[MAX_DIM];  // stride length for each index
  X                   mData[I*J*K*L*M*N]; // data container
  bool                mPeriodic=false;    // if true: disable bounds-check where possible

public:

  // return size without constructing
  static size_t Size();

  // constructor: allocate, don't initialize
  array();

  // constructor: copy from own class (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  array(const cppmat::tiny::array<U,RANK,I,J,K,L,M,N> &A);

  // constructor: copy from {...}
  array(const std::initializer_list<X> &A);

  // constructor: copy from dynamic size
  array(const cppmat::array<X> &A);

  // constructor: copy from view
  array(const cppmat::view::array<X,RANK,I,J,K,L,M,N> &A);

  // named constructor: initialize
  static array<X,RANK,I,J,K,L,M,N> Random  (X lower=(X)0, X upper=(X)1);
  static array<X,RANK,I,J,K,L,M,N> Arange  ();
  static array<X,RANK,I,J,K,L,M,N> Zero    ();
  static array<X,RANK,I,J,K,L,M,N> Ones    ();
  static array<X,RANK,I,J,K,L,M,N> Constant(X D);
  static array<X,RANK,I,J,K,L,M,N> Copy    (const std::vector<X> &D);

  // named constructor: copy
  template<typename It> static array<X,RANK,I,J,K,L,M,N> Copy(It first);
  template<typename It> static array<X,RANK,I,J,K,L,M,N> Copy(It first, It last);

  // return plain storage as vector
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  operator std::vector<U> () const;

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

  // bound check
  template<typename T> bool inBounds(T a) const;
  template<typename T> bool inBounds(T a, T b) const;
  template<typename T> bool inBounds(T a, T b, T c) const;
  template<typename T> bool inBounds(T a, T b, T c, T d) const;
  template<typename T> bool inBounds(T a, T b, T c, T d, T e) const;
  template<typename T> bool inBounds(T a, T b, T c, T d, T e, T f) const;

  // sign change
  array<X,RANK,I,J,K,L,M,N> operator- () const;
  array<X,RANK,I,J,K,L,M,N> operator+ () const;

  // arithmetic operators
  array<X,RANK,I,J,K,L,M,N>& operator*= (const array<X,RANK,I,J,K,L,M,N> &B);
  array<X,RANK,I,J,K,L,M,N>& operator/= (const array<X,RANK,I,J,K,L,M,N> &B);
  array<X,RANK,I,J,K,L,M,N>& operator+= (const array<X,RANK,I,J,K,L,M,N> &B);
  array<X,RANK,I,J,K,L,M,N>& operator-= (const array<X,RANK,I,J,K,L,M,N> &B);
  array<X,RANK,I,J,K,L,M,N>& operator*= (const       X                   &B);
  array<X,RANK,I,J,K,L,M,N>& operator/= (const       X                   &B);
  array<X,RANK,I,J,K,L,M,N>& operator+= (const       X                   &B);
  array<X,RANK,I,J,K,L,M,N>& operator-= (const       X                   &B);

  // absolute value
  array<X,RANK,I,J,K,L,M,N> abs() const;

  // norm (sum of absolute values)
  X norm() const;

  // return the indices that would sort an array
  array<size_t,RANK,I,J,K,L,M,N> argsort(bool ascending=true) const;

  // location of the minimum/maximum: plain storage (use decompress to convert to indices)
  size_t argmin() const;
  size_t argmax() const;

  // minimum
  X min() const;

  // maximum
  X max() const;

  // sum
  X sum() const;

  // mean
  double mean() const;

  // weighted average
  double average(const array<X,RANK,I,J,K,L,M,N> &weights, bool norm=true) const;

  // return array of booleans, based on condition
  array<int,RANK,I,J,K,L,M,N> equal        (const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> not_equal    (const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> greater      (const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> greater_equal(const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> less         (const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> less_equal   (const       X                   &D) const;
  array<int,RANK,I,J,K,L,M,N> equal        (const array<X,RANK,I,J,K,L,M,N> &D) const;
  array<int,RANK,I,J,K,L,M,N> not_equal    (const array<X,RANK,I,J,K,L,M,N> &D) const;
  array<int,RANK,I,J,K,L,M,N> greater      (const array<X,RANK,I,J,K,L,M,N> &D) const;
  array<int,RANK,I,J,K,L,M,N> greater_equal(const array<X,RANK,I,J,K,L,M,N> &D) const;
  array<int,RANK,I,J,K,L,M,N> less         (const array<X,RANK,I,J,K,L,M,N> &D) const;
  array<int,RANK,I,J,K,L,M,N> less_equal   (const array<X,RANK,I,J,K,L,M,N> &D) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

};

// =================================================================================================
// equality operators
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
bool operator!= (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
bool operator== (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// =================================================================================================
// external arithmetic operators
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator* (
  const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator/ (
  const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator+ (
  const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator- (
  const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator* (
  const array<X,RANK,I,J,K,L,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator/ (
  const array<X,RANK,I,J,K,L,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator+ (
  const array<X,RANK,I,J,K,L,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator- (
  const array<X,RANK,I,J,K,L,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator* (
  const X &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator/ (
  const X &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator+ (
  const X &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
array<X,RANK,I,J,K,L,M,N> operator- (
  const X &A, const array<X,RANK,I,J,K,L,M,N> &B);

// =================================================================================================
// print operator
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
std::ostream& operator<<(std::ostream& out, const array<X,RANK,I,J,K,L,M,N>& src);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

