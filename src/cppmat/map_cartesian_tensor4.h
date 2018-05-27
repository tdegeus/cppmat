/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_TENSOR4_H
#define CPPMAT_MAP_CARTESIAN_TENSOR4_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

// =================================================================================================
// cppmat::view::cartesian::tensor4
// =================================================================================================

template<class X, size_t ND>
class tensor4 : public cppmat::view::array<X,4,ND,ND,ND,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  tensor4();

  // constructor: map external pointer
  tensor4(const X *A);

  // get dimensions
  size_t ndim() const;

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

