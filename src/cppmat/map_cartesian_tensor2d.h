/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_TENSOR2D_H
#define CPPMAT_MAP_CARTESIAN_TENSOR2D_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

// =================================================================================================
// cppmat::view::cartesian::tensor2d
// =================================================================================================

template<typename X, size_t ND>
class tensor2d : public cppmat::view::diagonal::matrix<X,ND,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  tensor2d();

  // constructor: map external pointer
  tensor2d(const X *A);

  // named constructor: map external pointer
  static tensor2d<X,ND> Map(const X *D);

  // get dimensions
  size_t ndim() const;

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

