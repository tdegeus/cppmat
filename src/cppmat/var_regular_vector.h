/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_REGULAR_VECTOR_H
#define CPPMAT_VAR_REGULAR_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X>
class vector : public cppmat::array<X>
{
private:

  size_t N=0; // number of entries

private:

  // hide functions
  using cppmat::array<X>::reshape;
  using cppmat::array<X>::chrank;

public:

  // avoid name-hiding (also use the already defined overloads)
  using cppmat::array<X>::operator*=;
  using cppmat::array<X>::operator/=;
  using cppmat::array<X>::operator+=;
  using cppmat::array<X>::operator-=;

public:

  // constructor
  vector() = default;

  // constructor: allocate, don't initialize
  vector(size_t n);

  // constructor: copy
  vector(const array<X> &A);

  // constructor: copy
  vector(const std::vector<X> &D);

  // constructor: initialize
  static vector<X> Arange  (size_t n);
  static vector<X> Zero    (size_t n);
  static vector<X> Ones    (size_t n);
  static vector<X> Constant(size_t n, X D);

  // constructor: copy
  template<typename Itr> static vector<X> Copy(size_t n, Itr first);
  template<typename Itr> static vector<X> Copy(size_t n, Itr first, Itr last);

  // resize
  void resize(size_t n);

  // forward difference (x0, x1-x0, x2-x1, ...)
  vector<X> diff() const;

};

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

