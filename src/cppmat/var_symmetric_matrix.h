/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_SYMMETRIC_MATRIX_H
#define CPPMAT_VAR_SYMMETRIC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace symmetric {

// =================================================================================================
// cppmat::symmetric::matrix
// =================================================================================================

template<class X>
class matrix
{
protected:

  size_t              mSize=0;   // total size == data.size()
  static const size_t mRank=2;   // rank (number of axes)
  size_t              N=0;       // number of rows/columns
  std::vector<X>      mData;     // data container

public:

  // constructor: default
  matrix() = default;

  // constructor: allocate, don't initialize
  matrix(size_t m, size_t n);

  // constructor: copy
  matrix(const                   matrix<X> &A);
  matrix(const cppmat::diagonal::matrix<X> &A);

  // named constructor: initialize
  static matrix<X> Random  (size_t m, size_t n, X lower=(X)0, X upper=(X)1);
  static matrix<X> Arange  (size_t m, size_t n);
  static matrix<X> Zero    (size_t m, size_t n);
  static matrix<X> Ones    (size_t m, size_t n);
  static matrix<X> Constant(size_t m, size_t n, X D);
  static matrix<X> Copy    (size_t m, size_t n, const std::vector<X> &D);

  // named constructor: copy
  template<typename Itr> static matrix<X> Copy     (size_t m, size_t n, Itr first);
  template<typename Itr> static matrix<X> Copy     (size_t m, size_t n, Itr first, Itr last);
  template<typename Itr> static matrix<X> CopyDense(size_t m, size_t n, Itr first);
  template<typename Itr> static matrix<X> CopyDense(size_t m, size_t n, Itr first, Itr last);

  // return plain storage as vector
  std::vector<X> asVector() const;

  // resize
  void resize(size_t m, size_t n);

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
  X&       operator()(size_t a, size_t b);
  const X& operator()(size_t a, size_t b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

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
  auto item(size_t a, size_t b);
  auto item(size_t a, size_t b) const;

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
  matrix<X> operator- () const;
  matrix<X> operator+ () const;

  // arithmetic operators
  matrix<X>& operator*= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator/= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator+= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator-= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator*= (const                           X  &B);
  matrix<X>& operator/= (const                           X  &B);
  matrix<X>& operator+= (const                           X  &B);
  matrix<X>& operator-= (const                           X  &B);

  // extra arithmetic operators
  matrix<X>& operator*= (const cppmat::diagonal ::matrix<X> &B);
  matrix<X>& operator+= (const cppmat::diagonal ::matrix<X> &B);
  matrix<X>& operator-= (const cppmat::diagonal ::matrix<X> &B);

  // absolute value
  matrix<X> abs() const;

  // norm (sum of absolute values)
  X norm() const;

  // return the indices that would sort the matrix
  matrix<size_t> argsort(bool ascending=true) const;

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
  double average(const matrix<X> &weights, bool norm=true) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;

  // find the plain storage indices of all entries equal to some constant
  std::vector<size_t> where(X D) const;

};

// external arithmetic operators (cppmat::symmetric::matrix)
template<class X> inline matrix<X> operator* (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator- (const matrix<X> &A, const matrix<X> &B);
template<class X> inline matrix<X> operator* (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator/ (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator+ (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator- (const matrix<X> &A, const        X  &B);
template<class X> inline matrix<X> operator* (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator/ (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator+ (const        X  &A, const matrix<X> &B);
template<class X> inline matrix<X> operator- (const        X  &A, const matrix<X> &B);

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

