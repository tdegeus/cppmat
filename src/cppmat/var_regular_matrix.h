/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_MATRIX_H
#define CPPMAT_VAR_REGULAR_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template<class X>
class matrix : public cppmat::array<X>
{
protected:

  size_t M=0; // number of rows
  size_t N=0; // number of columns

private:

  // hide functions
  using cppmat::array<X>::chrank;

public:

  // avoid name-hiding (also use the already defined overloads)
  using cppmat::array<X>::operator*=;
  using cppmat::array<X>::operator/=;
  using cppmat::array<X>::operator+=;
  using cppmat::array<X>::operator-=;

public:

  // constructor
  matrix() = default;

  // constructor: allocate, don't initialize
  matrix(size_t m, size_t n);

  // constructor: copy
  matrix(const cppmat::array<X> &A);

  // constructor: copy
  #ifndef CPPMAT_NOCONVERT
  matrix(const cppmat::symmetric::matrix<X> &A);
  matrix(const cppmat::diagonal ::matrix<X> &A);
  #endif

  // constructor: copy
  matrix(size_t m, size_t n, const std::vector<X> &D);

  // constructor: initialize
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

  // iterator to the first and last entry of a row
  auto beginRow(size_t i);
  auto beginRow(size_t i) const;
  auto endRow  (size_t i);
  auto endRow  (size_t i) const;

  // extra arithmetic operators
  matrix<X>& operator*= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator/= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator+= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator-= (const cppmat::symmetric::matrix<X> &B);
  matrix<X>& operator*= (const cppmat::diagonal ::matrix<X> &B);
  matrix<X>& operator+= (const cppmat::diagonal ::matrix<X> &B);
  matrix<X>& operator-= (const cppmat::diagonal ::matrix<X> &B);
};

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

