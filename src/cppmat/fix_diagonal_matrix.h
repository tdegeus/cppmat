/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_DIAGONAL_MATRIX_H
#define CPPMAT_FIX_DIAGONAL_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace diagonal {

// =================================================================================================
// cppmat::tiny::diagonal::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix
{
  static_assert( N == M , "Must be square" );

protected:

  static const size_t mSize=N;          // total size == data.size()
  static const size_t mRank=2;          // rank (number of axes)
  X                   mData[N];         // data container
  X                   mZero[1];         // pointer to a zero entry
  bool                mPeriodic=false;  // if true: disable bounds-check where possible

public:

  // return size without constructing
  static size_t Size();

  // constructor: allocate, don't initialize
  matrix();

  // constructor: copy from own class (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  matrix(const cppmat::tiny::diagonal::matrix<U,M,N> &A);

  // constructor: copy from dynamic size
  matrix(const cppmat::diagonal::matrix<X> &A);

  // constructor: copy from view
  matrix(const cppmat::view::diagonal::matrix<X,M,N> &A);

  // named constructor: initialize
  static matrix<X,M,N> Random  (X lower=(X)0, X upper=(X)1);
  static matrix<X,M,N> Arange  ();
  static matrix<X,M,N> Zero    ();
  static matrix<X,M,N> Ones    ();
  static matrix<X,M,N> Constant(X D);
  static matrix<X,M,N> Copy    (const std::vector<X> &D);

  // named constructor: copy
  template<typename Itr> static matrix<X,M,N> Copy     (Itr first);
  template<typename Itr> static matrix<X,M,N> Copy     (Itr first, Itr last);
  template<typename Itr> static matrix<X,M,N> CopyDense(Itr first);
  template<typename Itr> static matrix<X,M,N> CopyDense(Itr first, Itr last);

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

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  X&       operator()(int a, int b);
  const X& operator()(int a, int b) const;

  // index operators: access using matrix-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  X&       operator()(T a, T b);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(int a, int b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
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

  // iterator to specific entry: access using matrix-indices
  auto item(int a, int b);
  auto item(int a, int b) const;

  // iterator to specific entry: access using matrix-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b);

  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b) const;

  // initialization
  void setRandom(X lower=(X)0, X upper=(X)1);
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy     (Iterator first);
  template<typename Iterator> void setCopy     (Iterator first, Iterator last);
  template<typename Iterator> void setCopyDense(Iterator first);
  template<typename Iterator> void setCopyDense(Iterator first, Iterator last);

  // copy to target
  template<typename Iterator> void copyTo     (Iterator first) const;
  template<typename Iterator> void copyTo     (Iterator first, Iterator last) const;
  template<typename Iterator> void copyToDense(Iterator first) const;
  template<typename Iterator> void copyToDense(Iterator first, Iterator last) const;

  // sign change
  matrix<X,M,N> operator- () const;
  matrix<X,M,N> operator+ () const;

  // arithmetic operators
  matrix<X,M,N>& operator*= (const cppmat::tiny::diagonal::matrix<X,M,N> &B);
  matrix<X,M,N>& operator+= (const cppmat::tiny::diagonal::matrix<X,M,N> &B);
  matrix<X,M,N>& operator-= (const cppmat::tiny::diagonal::matrix<X,M,N> &B);
  matrix<X,M,N>& operator*= (const                                X      &B);
  matrix<X,M,N>& operator/= (const                                X      &B);

  // extra arithmetic operators
  matrix<X,M,N>& operator*= (const cppmat::tiny::           matrix<X,M,N> &B);
  matrix<X,M,N>& operator/= (const cppmat::tiny::           matrix<X,M,N> &B);
  matrix<X,M,N>& operator*= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);
  matrix<X,M,N>& operator/= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);

  // absolute value
  matrix<X,M,N> abs() const;

  // norm (sum of absolute values)
  X norm() const;

  // return the indices that would sort the matrix
  matrix<size_t,M,N> argsort(bool ascending=true) const;

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
  double average(const matrix<X,M,N> &weights, bool norm=true) const;

  // return matrix of booleans, based on condition
  matrix<int,M,N> equal        (const        X      &D) const;
  matrix<int,M,N> not_equal    (const        X      &D) const;
  matrix<int,M,N> greater      (const        X      &D) const;
  matrix<int,M,N> greater_equal(const        X      &D) const;
  matrix<int,M,N> less         (const        X      &D) const;
  matrix<int,M,N> less_equal   (const        X      &D) const;
  matrix<int,M,N> equal        (const matrix<X,M,N> &D) const;
  matrix<int,M,N> not_equal    (const matrix<X,M,N> &D) const;
  matrix<int,M,N> greater      (const matrix<X,M,N> &D) const;
  matrix<int,M,N> greater_equal(const matrix<X,M,N> &D) const;
  matrix<int,M,N> less         (const matrix<X,M,N> &D) const;
  matrix<int,M,N> less_equal   (const matrix<X,M,N> &D) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

};

// =================================================================================================
// external arithmetic operators (cppmat::tiny::diagonal)
// =================================================================================================

template<class X, size_t M, size_t N>
matrix<X,M,N> operator* (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
matrix<X,M,N> operator+ (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
matrix<X,M,N> operator- (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
matrix<X,M,N> operator* (const matrix<X,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
matrix<X,M,N> operator/ (const matrix<X,M,N> &A, const X &B);

// -------------------------------------------------------------------------------------------------

template<class X, size_t M, size_t N>
matrix<X,M,N> operator* (const X &A, const matrix<X,M,N> &B);

// =================================================================================================
// print operator
// =================================================================================================

template<class X, size_t M, size_t N>
std::ostream& operator<<(std::ostream& out, const matrix<X,M,N>& src);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

