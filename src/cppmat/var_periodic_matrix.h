/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_PERIODIC_MATRIX_H
#define CPPMAT_VAR_PERIODIC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template<class X>
class matrix : public cppmat::periodic::array<X>
{
protected:

  int M=0; // number of rows
  int N=0; // number of columns

private:

  // hide functions
  using cppmat::periodic::array<X>::chrank;

public:

  // constructor
  matrix() = default;

  // constructor: allocate, don't initialize
  matrix(size_t m, size_t n);

  // constructor: copy
  matrix(const cppmat::periodic::array<X> &A);

  // constructor: copy
  matrix(const cppmat::array<X> &A);

  // constructor: copy
  matrix(size_t m, size_t n, const std::vector<X> &D);

  // constructor: initialize
  static matrix<X> Random  (size_t m, size_t n, X lower=(X)0, X upper=(X)1);
  static matrix<X> Arange  (size_t m, size_t n);
  static matrix<X> Zero    (size_t m, size_t n);
  static matrix<X> Ones    (size_t m, size_t n);
  static matrix<X> Constant(size_t m, size_t n, X D);

  // constructor: copy
  template<typename Itr> static matrix<X> Copy(size_t m, size_t n, Itr first);
  template<typename Itr> static matrix<X> Copy(size_t m, size_t n, Itr first, Itr last);

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

  // get dimensions
  size_t rows() const;
  size_t cols() const;

  // iterator to the first and last entry of a row
  auto beginRow(int i);
  auto beginRow(int i) const;
  auto endRow  (int i);
  auto endRow  (int i) const;

  // iterator to the first and last entry of a row
  auto beginRow(size_t i);
  auto beginRow(size_t i) const;
  auto endRow  (size_t i);
  auto endRow  (size_t i) const;

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

