/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_H
#define CPPMAT_VAR_CARTESIAN_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor4<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor4<X> &B
);

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2<X> &B
);

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2s<X> &B
);

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2d<X> &B
);

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor4<X> T(const cppmat::cartesian::tensor4<X> &A);

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> RT(const cppmat::cartesian::tensor4<X> &A);

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> LT(const cppmat::cartesian::tensor4<X> &A);

// =================================================================================================

}} // namespace ...

#endif

