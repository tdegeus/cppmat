/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_REGULAR_MATRIX_H
#define CPPMAT_MAP_REGULAR_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

// =================================================================================================
// cppmat::view::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix : public cppmat::view::array<X,2,M,N>
{
public:

  // constructor: allocate, don't initialize
  matrix();

  // constructor: map external pointer
  matrix(const X *A);

  // named constructor: map external pointer
  static matrix<X,M,N> Map(const X *D);

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

