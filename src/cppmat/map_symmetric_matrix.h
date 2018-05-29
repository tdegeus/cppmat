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

template<class X, size_t M, size_t N>
class matrix
{
  static_assert( N == M , "Must be square" );

protected:

  static const size_t mSize=(N+1)*N/2;  // total size == data.size()
  static const size_t mRank=2;          // rank (number of axes)
  const X            *mData;            // data container

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
  std::vector<X> asVector() const;

  // get dimensions
  size_t size() const;
  size_t rank() const;
  size_t shape(int    i) const;
  size_t shape(size_t i) const;
  std::vector<size_t> shape() const;

  // index operators: access plain storage
  const X& operator[](size_t i) const;

  // index operators: access using matrix-indices
  const X& operator()(size_t a, size_t b) const;

  // index operators: matrix-indices -> plain storage (a,b -> i)
  size_t compress(size_t a, size_t b) const;

  // index operators: plain storage -> matrix-indices (i -> a,b)
  std::vector<size_t> decompress(size_t i) const;

  // pointer to data
  const X* data() const;

  // iterator to first and last entry
  auto begin() const;
  auto end() const;

  // iterator to specific entry: access plain storage
  auto index(size_t i) const;

  // iterator to specific entry: access using matrix-indices
  auto item(size_t a, size_t b) const;

  // initialization
  void setMap(const X *D);

  // copy to target
  template<typename Iterator> void copyTo     (Iterator first) const;
  template<typename Iterator> void copyTo     (Iterator first, Iterator last) const;
  template<typename Iterator> void copyToDense(Iterator first) const;
  template<typename Iterator> void copyToDense(Iterator first, Iterator last) const;

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

  // find the plain storage indices of all entries equal to some constant
  std::vector<size_t> where(X D) const;

};

// print operator
template<class X, size_t M, size_t N>
std::ostream& operator<<(std::ostream& out, const matrix<X,M,N>& src);

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
