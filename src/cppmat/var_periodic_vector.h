/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_PERIODIC_VECTOR_H
#define CPPMAT_VAR_PERIODIC_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::vector
// =================================================================================================

template<class X>
class vector : public cppmat::periodic::array<X>
{
private:

  // hide functions
  using cppmat::periodic::array<X>::reshape;
  using cppmat::periodic::array<X>::chrank;

public:

  // constructor: default
  vector() = default;

  // constructor: allocate, don't initialize
  vector(size_t n);

  // constructor: copy
  vector(const cppmat::array <X> &A);
  vector(const std   ::vector<X> &A);

  // constructor: copy
  template<size_t n> vector(const cppmat::tiny::periodic::vector<X,n> &A);
  template<size_t n> vector(const cppmat::view::periodic::vector<X,n> &A);

  // named constructor: initialize
  static vector<X> Random  (size_t n, X lower=(X)0, X upper=(X)1);
  static vector<X> Arange  (size_t n);
  static vector<X> Zero    (size_t n);
  static vector<X> Ones    (size_t n);
  static vector<X> Constant(size_t n, X D);
  static vector<X> Copy    (size_t n, const std::vector<X> &D);
  static vector<X> Copy    (          const std::vector<X> &D);

  // named constructor: copy
  template<typename Itr> static vector<X> Copy(size_t n, Itr first);
  template<typename Itr> static vector<X> Copy(size_t n, Itr first, Itr last);

  // resize
  void resize(size_t n);

  // forward difference (x0, x1-x0, x2-x1, ...)
  vector<X> diff() const;

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

