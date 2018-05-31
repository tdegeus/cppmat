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

  // constructor: default
  matrix() = default;

  // constructor: allocate, don't initialize
  matrix(size_t m, size_t n);

  // constructor: copy from parent (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  matrix(const cppmat::array<U> &A);

  // constructor: copy from other class
  matrix(const cppmat::symmetric::matrix<X> &A);
  matrix(const cppmat::diagonal ::matrix<X> &A);

  // constructor: copy from fixed size
  template<size_t m, size_t n> matrix(const cppmat::tiny::matrix<X,m,n> &A);

  // constructor: copy from view
  template<size_t m, size_t n> matrix(const cppmat::view::matrix<X,m,n> &A);

  // named constructor: initialize
  static matrix<X> Random  (size_t m, size_t n, X lower=(X)0, X upper=(X)1);
  static matrix<X> Arange  (size_t m, size_t n);
  static matrix<X> Zero    (size_t m, size_t n);
  static matrix<X> Ones    (size_t m, size_t n);
  static matrix<X> Constant(size_t m, size_t n, X D);
  static matrix<X> Copy    (size_t m, size_t n, const std::vector<X> &D);

  // named constructor: copy
  template<typename Itr> static matrix<X> Copy(size_t m, size_t n, Itr first);
  template<typename Itr> static matrix<X> Copy(size_t m, size_t n, Itr first, Itr last);

  // resize
  void resize (size_t m, size_t n);
  void reshape(size_t m, size_t n);

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

