/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_REGULAR_ARRAY_H
#define CPPMAT_MAP_REGULAR_ARRAY_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

// =================================================================================================
// cppmat::view::array
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
class array
{
protected:

  static const size_t MAX_DIM=6;          // maximum number of dimensions
  static const size_t mSize=I*J*K*L*M*N;  // total size == data.size() == prod(shape)
  static const size_t mRank=RANK;         // rank (number of axes)
  size_t              mShape  [MAX_DIM];  // number of entries along each axis
  size_t              mStrides[MAX_DIM];  // stride length for each index
  const X            *mData;              // data container
  bool                mPeriodic=false;    // if true: disable bounds-check where possible

public:

  // return size without constructing
  static size_t Size();

  // constructor: allocate, don't initialize
  array();

  // constructor: map external pointer
  array(const X *A);

  // named constructor: map external pointer
  static array<X,RANK,I,J,K,L,M,N> Map(const X *D);

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

  // get using a different return type
  template<typename U> U size() const;
  template<typename U> U rank() const;
  template<typename U> U shape(int    i) const;
  template<typename U> U shape(size_t i) const;
  template<typename U> std::vector<U> shape() const;
  template<typename U> std::vector<U> strides(bool bytes=false) const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using array-indices
  const X& operator()(int a) const;
  const X& operator()(int a, int b) const;
  const X& operator()(int a, int b, int c) const;
  const X& operator()(int a, int b, int c, int d) const;
  const X& operator()(int a, int b, int c, int d, int e) const;
  const X& operator()(int a, int b, int c, int d, int e, int f) const;

  // index operators: access using array-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d, T e) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b, T c, T d, T e, T f) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of array-indices (a,b,c,...)
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
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using array-indices
  auto item(int a) const;
  auto item(int a, int b) const;
  auto item(int a, int b, int c) const;
  auto item(int a, int b, int c, int d) const;
  auto item(int a, int b, int c, int d, int e) const;
  auto item(int a, int b, int c, int d, int e, int f) const;

  // iterator to specific entry: access using array-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e) const;

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b, T c, T d, T e, T f) const;

  // initialization
  void setMap(const X *D);

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

  // norm (sum of absolute values)
  X norm() const;

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

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;
  size_t              where(int    index) const;
  size_t              where(size_t index) const;

};

// =================================================================================================
// equality operators
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
bool operator!= (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
bool operator== (const array<X,RANK,I,J,K,L,M,N> &A, const array<X,RANK,I,J,K,L,M,N> &B);

// =================================================================================================
// print operator
// =================================================================================================

template<typename X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
std::ostream& operator<<(std::ostream& out, const array<X,RANK,I,J,K,L,M,N>& src);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

