/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VECTOR_H
#define CPPMAT_VECTOR_H

#include "cppmat.h"

namespace cppmat {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X>
class vector
{
private:

  std::vector<X> m_data;    // data container
  size_t         m_n=0;     // number of columns
  size_t         m_size=0;  // total size

public:

  // constructors
  vector(){};
  vector(size_t n);
  vector(size_t n, X D);

  template<typename Iterator>
  vector(Iterator first, Iterator last);

  // resize
  void resize(size_t n);

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
  vector<X>& operator*= (const vector<X> &B);
  vector<X>& operator/= (const vector<X> &B);
  vector<X>& operator+= (const vector<X> &B);
  vector<X>& operator-= (const vector<X> &B);
  vector<X>& operator*= (const        X  &B);
  vector<X>& operator/= (const        X  &B);
  vector<X>& operator+= (const        X  &B);
  vector<X>& operator-= (const        X  &B);

  // basic algebra
  X      min() const;
  X      max() const;
  X      sum() const;
  double mean() const;
  double average(const vector<X> &weights) const;

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

} // namespace ...

#endif

