/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MAP_CARTESIAN_TENSOR2S_H
#define CPPMAT_MAP_CARTESIAN_TENSOR2S_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

// =================================================================================================
// cppmat::view::cartesian::tensor2s
// =================================================================================================

template<typename X, size_t ND>
class tensor2s : public cppmat::view::symmetric::matrix<X,ND,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  tensor2s();

  // constructor: map external pointer
  tensor2s(const X *A);

  // named constructor: map external pointer
  static tensor2s<X,ND> Map(const X *D);

  // get dimensions
  size_t ndim() const;

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

