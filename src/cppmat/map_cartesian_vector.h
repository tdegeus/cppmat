/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_VECTOR_H
#define CPPMAT_MAP_CARTESIAN_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

// =================================================================================================
// cppmat::view::cartesian::vector
// =================================================================================================

template<typename X, size_t ND>
class vector : public cppmat::view::vector<X,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  vector();

  // constructor: map external pointer
  vector(const X *A);

  // named constructor: map external pointer
  static vector<X,ND> Map(const X *D);

  // get dimensions
  size_t ndim() const;

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

