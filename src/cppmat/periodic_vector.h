/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_VECTOR_H
#define CPPMAT_PERIODIC_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X>
class vector
{
private:

  std::vector<X> mData;    // data container
  size_t         mSize=0;  // total size
  size_t         N=0;      // number of columns
  int            NI=0;     // == N, but int

public:

  // constructor
  vector() = default;
  vector(size_t n);

  // constructor: initialize
  static vector<X> Arange  (size_t n);
  static vector<X> Zero    (size_t n);
  static vector<X> Ones    (size_t n);
  static vector<X> Constant(size_t n, X D);

  // constructor: initialize by copying from external object
  template<typename Iterator>
  static vector<X> Copy(Iterator first, Iterator last);

  // resize
  void resize(size_t n);

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
  X&       operator()(int a);
  const X& operator()(int a) const;

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
  auto item(int a);
  auto item(int a) const;

  // basic initialization
  void setArange();
  void setZero();
  void setOnes();
  void setConstant(X D);
  template<typename Iterator> void setCopy(Iterator first, Iterator last);

  // arithmetic operators
  vector<X>& operator*= (const vector<X> &B);
  vector<X>& operator/= (const vector<X> &B);
  vector<X>& operator+= (const vector<X> &B);
  vector<X>& operator-= (const vector<X> &B);
  vector<X>& operator*= (const        X  &B);
  vector<X>& operator/= (const        X  &B);
  vector<X>& operator+= (const        X  &B);
  vector<X>& operator-= (const        X  &B);

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
  double average(const vector<X> &weights, bool norm=true) const;

  // find all non-zero entries
  vector<size_t> where() const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X> inline vector<X> operator* (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator/ (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator+ (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator- (const vector<X> &A, const vector<X> &B);
template<class X> inline vector<X> operator* (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator/ (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator+ (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator- (const vector<X> &A, const        X  &B);
template<class X> inline vector<X> operator* (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator/ (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator+ (const        X  &A, const vector<X> &B);
template<class X> inline vector<X> operator- (const        X  &A, const vector<X> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

