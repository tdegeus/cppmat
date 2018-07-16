/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_REGULAR_VECTOR_H
#define CPPMAT_MAP_REGULAR_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

// =================================================================================================
// cppmat::view::vector
// =================================================================================================

template<typename X, size_t N>
class vector : public cppmat::view::array<X,1,N>
{
public:

  // constructor: allocate, don't initialize
  vector();

  // constructor: map external pointer
  vector(const X *A);

  // named constructor: map external pointer
  static vector<X,N> Map(const X *D);

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

