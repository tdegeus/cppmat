/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_SYMMETRIC_MATRIX_H
#define CPPMAT_MAP_SYMMETRIC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace symmetric {

// =================================================================================================
// cppmat::view::symmetric::matrix
// =================================================================================================

template<typename X, size_t M, size_t N>
class matrix
{
  static_assert( N == M , "Must be square" );

protected:

  static const size_t mSize=(N+1)*N/2;  // total size == data.size()
  static const size_t mRank=2;          // rank (number of axes)
  const X            *mData;            // data container
  bool                mPeriodic=false;  // if true: disable bounds-check where possible

public:

  // return size without constructing
  static size_t Size();

  // constructor: allocate, don't initialize
  matrix();

  // constructor: map external pointer
  matrix(const X *A);

  // named constructor: map external pointer
  static matrix<X,M,N> Map(const X *D);

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

  // get using a different return type
  template<typename U> U size() const;
  template<typename U> U rank() const;
  template<typename U> U shape(int    i) const;
  template<typename U> U shape(size_t i) const;
  template<typename U> std::vector<U> shape() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  const X& operator()(int a, int b) const;

  // index operators: access using matrix-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  const X& operator()(T a, T b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(int a, int b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  size_t compress(T a, T b) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
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

  // iterator to specific entry: access using matrix-indices
  auto item(int a, int b) const;

  // iterator to specific entry: access using matrix-indices
  template<typename T, typename=typename std::enable_if<std::is_unsigned<T>::value,void>::type>
  auto item(T a, T b) const;

  // initialization
  void setMap(const X *D);

  // copy to target
  template<typename Iterator> void copyTo     (Iterator first) const;
  template<typename Iterator> void copyTo     (Iterator first, Iterator last) const;
  template<typename Iterator> void copyToDense(Iterator first) const;
  template<typename Iterator> void copyToDense(Iterator first, Iterator last) const;

  // check in an index in within the matrix bound:
  // ( a >= 0 and a < shape[0] ) or ( periodic )
  template<typename T> bool inBounds(T a) const;
  template<typename T> bool inBounds(T a, T b) const;

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
  double average(const matrix<X,M,N> &weights, bool norm=true) const;

  // find the plain storage indices of all non-zero entries
  std::vector<size_t> where() const;
  size_t              where(int    index) const;
  size_t              where(size_t index) const;

};

// =================================================================================================
// equality operators
// =================================================================================================

template<typename X, size_t M, size_t N>
bool operator!= (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

// -------------------------------------------------------------------------------------------------

template<typename X, size_t M, size_t N>
bool operator== (const matrix<X,M,N> &A, const matrix<X,M,N> &B);

// =================================================================================================
// print operator
// =================================================================================================

template<typename X, size_t M, size_t N>
std::ostream& operator<<(std::ostream& out, const matrix<X,M,N>& src);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

