/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_PERIODIC_MATRIX_H
#define CPPMAT_FIX_PERIODIC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace periodic {

// =================================================================================================
// cppmat::tiny::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix : public cppmat::tiny::periodic::array<X,2,M,N>
{
public:

  // constructor: allocate, don't initialize
  matrix();

  // constructor: copy from parent
  matrix(const cppmat::tiny::array<X,2,M,N> &A);

  // constructor: copy from dynamic size
  matrix(const cppmat::periodic::matrix<X> &A);

  // constructor: copy from view
  matrix(const cppmat::view::periodic::matrix<X,M,N> &A);

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

