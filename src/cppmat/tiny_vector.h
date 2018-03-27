/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TINY_VECTOR_H
#define CPPMAT_TINY_VECTOR_H

#include "cppmat.h"

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::vector
// =================================================================================================

template<class X, size_t n>
class vector
{
private:

  // data container
  typename std::remove_const<X>::type m_container[n];
  // pointer to container (may point outside)
  X *m_data;
  // total number of entries
  size_t m_size=n;

public:

  // constructors
  vector();
  vector(X D);

  template<typename Iterator>
  vector(Iterator first, Iterator last);

  // copy constructor
  vector(const vector<X,n> &D);

  // assignment operator
  vector<X,n>& operator= (const vector<X,n> &D);

  // map external pointer
  void map(X *D);

  // get dimensions
  size_t size() const;
  size_t ndim() const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;
  std::vector<size_t> strides(bool bytes=false) const;

  // index operators: access plain storage
  X&       operator[](size_t i);
  const X& operator[](size_t i) const;

  // index operators: access using matrix indices
  X&       operator()(size_t a);
  const X& operator()(size_t a) const;

  // pointer / iterators
  X*       data();
  const X* data() const;
  auto     begin();
  auto     begin() const;
  auto     end();
  auto     end() const;

  // basic initialization
  void setConstant(X D);
  void setZero();
  void setOnes();
  void zeros();
  void ones();

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
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const vector<X,n> &weights) const;

  // formatted print; NB also "operator<<" is defined
  void printf(std::string fmt) const;

};

// arithmetic operators
template<class X, size_t n> inline vector<X,n> operator* (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator/ (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator+ (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator- (const vector<X,n> &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator* (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator/ (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator+ (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator- (const vector<X,n> &A, const        X    &B);
template<class X, size_t n> inline vector<X,n> operator* (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator/ (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator+ (const        X    &A, const vector<X,n> &B);
template<class X, size_t n> inline vector<X,n> operator- (const        X    &A, const vector<X,n> &B);



// =================================================================================================

}} // namespace ...

#endif

