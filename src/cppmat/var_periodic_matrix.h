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
private:

  // hide functions
  using cppmat::periodic::array<X>::chrank;

public:

  // constructor: default
  matrix() = default;

  // constructor: allocate, don't initialize
  matrix(size_t m, size_t n);

  // constructor: copy
  matrix(const cppmat::array<X> &A);

  // constructor: copy
  template<size_t m, size_t n> matrix(const cppmat::tiny::periodic::matrix<X,m,n> &A);
  template<size_t m, size_t n> matrix(const cppmat::view::periodic::matrix<X,m,n> &A);

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

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

