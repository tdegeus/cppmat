/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_VECTOR_H
#define CPPMAT_TINY_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::vector
// =================================================================================================

template<class X, size_t n>
class vector
{
private:

  static const size_t m_size=n;  // total size
  static const size_t m_n=n;     // number of columns
  X m_data[m_size];              // data container

public:

  // constructor
  vector(){};

  // constructor: initialize
  static vector<X,n> Arange();
  static vector<X,n> Zero();
  static vector<X,n> Ones();
  static vector<X,n> Constant(X D);

  // constructor: initialize by copying from external object
  template<typename Iterator> static vector<X,n> Copy(Iterator first);
  template<typename Iterator> static vector<X,n> Copy(Iterator first, Iterator last);

  // information without constructing
  static size_t Size();

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

  // index operators: access using list-index
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;

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

  // iterator to specific entry: access using list-index
  auto item(size_t a);
  auto item(size_t a) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

  // arithmetic operators
  vector<X,n>& operator*= (const vector<X,n> &B);
  vector<X,n>& operator/= (const vector<X,n> &B);
  vector<X,n>& operator+= (const vector<X,n> &B);
  vector<X,n>& operator-= (const vector<X,n> &B);
  vector<X,n>& operator*= (const        X    &B);
  vector<X,n>& operator/= (const        X    &B);
  vector<X,n>& operator+= (const        X    &B);
  vector<X,n>& operator-= (const        X    &B);

  // basic algebra
  // - absolute value
  void abs();
  // - location of the minimum/maximum
  size_t argmin() const;
  size_t argmax() const;
  // - minimum
  X minCoeff() const;
  // - maximum
  X maxCoeff() const;
  // - sum
  X sum() const;
  // - mean
  double mean() const;
  // - weighted average
  double average(const vector<X,n> &weights, bool norm=true) const;

  // find all non-zero entries
  cppmat::vector<size_t> where() const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator* (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline vector<X,n> operator/ (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline vector<X,n> operator+ (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline vector<X,n> operator- (const vector<X,n> &A, const        X    &B);

template<class X, size_t n>
inline vector<X,n> operator* (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator/ (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator+ (const        X    &A, const vector<X,n> &B);

template<class X, size_t n>
inline vector<X,n> operator- (const        X    &A, const vector<X,n> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

