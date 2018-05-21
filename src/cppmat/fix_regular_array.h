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
// cppmat::array
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J=1, size_t K=1, size_t L=1, size_t M=1, size_t N=1>
class array
{
protected:

  static const size_t MAX_DIM=6;  // maximum number of dimensions
  size_t mSize=I*J*K*L*M*N;       // total size == data.size() == prod(shape)
  size_t mRank=RANK;              // rank (number of axes)
  size_t mShape  [MAX_DIM];       // number of entries along each axis
  size_t mStrides[MAX_DIM];       // stride length for each index
  X      mData[I*J*K*L*M*N];      // data container

public:

  // return size without constructing
  static size_t Size();

  // constructor: allocate, don't initialize
  array();

  // constructor: copy
  array(const array<X,RANK,I,J,K,L,M,N> &A);
  array(const cppmat::array<X>          &A);

  // named constructor: initialize
  static array<X,RANK,I,J,K,L,M,N> Random  (X lower=(X)0, X upper=(X)1);
  static array<X,RANK,I,J,K,L,M,N> Arange  ();
  static array<X,RANK,I,J,K,L,M,N> Zero    ();
  static array<X,RANK,I,J,K,L,M,N> Ones    ();
  static array<X,RANK,I,J,K,L,M,N> Constant(X D);

  // named constructor: copy
  static array<X,RANK,I,J,K,L,M,N> Copy(const std::vector<X> &D);

  // named constructor: copy
  template<typename It> static array<X,RANK,I,J,K,L,M,N> Copy(It first);
  template<typename It> static array<X,RANK,I,J,K,L,M,N> Copy(It first, It last);

  // copy constructor
  operator cppmat::array<X> () const;

  // return plain storage as vector
  std::vector<X> asVector() const;

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

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

  // find the plain storage indices of all entries equal to some constant
  std::vector<size_t> where(X D) const;

};

// external arithmetic operators
template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (array<X,RANK,I,J,K,L,M,N> A, const array<X,RANK,I,J,K,L,M,N> &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (array<X,RANK,I,J,K,L,M,N> A, const X &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (array<X,RANK,I,J,K,L,M,N> A, const X &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (array<X,RANK,I,J,K,L,M,N> A, const X &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (array<X,RANK,I,J,K,L,M,N> A, const X &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator* (const X &A,       array<X,RANK,I,J,K,L,M,N>  B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator/ (const X &A, const array<X,RANK,I,J,K,L,M,N> &B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator+ (const X &A,       array<X,RANK,I,J,K,L,M,N>  B);

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
inline
array<X,RANK,I,J,K,L,M,N> operator- (const X &A, const array<X,RANK,I,J,K,L,M,N> &B);


// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

