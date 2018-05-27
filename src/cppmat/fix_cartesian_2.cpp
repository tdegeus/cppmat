/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_2_CPP
#define CPPMAT_FIX_CARTESIAN_2_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// tensor products: ddot
// =================================================================================================

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2<X,2> &A, const cppmat::tiny::cartesian::tensor2d<X,2> &B
)
{
  return A[0]*B[0] + A[3]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,2> &A, const cppmat::tiny::cartesian::tensor2s<X,2> &B
)
{
  X C;

  C  = A[0] * B[0];
  C += A[1] * B[1] * static_cast<X>(2);
  C += A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,2> &A, const cppmat::tiny::cartesian::tensor2d<X,2> &B
)
{
  return A[0]*B[0] + A[2]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,2> &A, const cppmat::tiny::cartesian::tensor2<X,2> &B
)
{
  return A[0]*B[0] + A[1]*B[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,2> &A, const cppmat::tiny::cartesian::tensor2s<X,2> &B
)
{
  return A[0]*B[0] + A[1]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,2> &A, const cppmat::tiny::cartesian::tensor2d<X,2> &B
)
{
  return A[0]*B[0] + A[1]*B[1];
}

// =================================================================================================
// determinant
// =================================================================================================

template<class X>
CPPMAT_INLINE
X det(const cppmat::tiny::cartesian::tensor2<X,2> &A)
{
 return A[0] * A[3] - A[1] * A[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
X det(const cppmat::tiny::cartesian::tensor2s<X,2> &A)
{
 return A[0] * A[2] - A[1] * A[1];
}

// =================================================================================================
// inverse
// =================================================================================================

template<class X>
CPPMAT_INLINE
cppmat::tiny::cartesian::tensor2<X,2> inv(const cppmat::tiny::cartesian::tensor2<X,2> &A)
{
  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::tiny::cartesian::tensor2<X,2> C;

  // compute inverse
  C[0] =                      A[3] / det;
  C[1] = static_cast<X>(-1) * A[1] / det;
  C[2] = static_cast<X>(-1) * A[2] / det;
  C[3] =                      A[0] / det;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
CPPMAT_INLINE
cppmat::tiny::cartesian::tensor2s<X,2> inv(const cppmat::tiny::cartesian::tensor2s<X,2> &A)
{
  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::tiny::cartesian::tensor2s<X,2> C;

  // compute inverse
  C[0] =                      A[2] / det;
  C[1] = static_cast<X>(-1) * A[1] / det;
  C[2] =                      A[0] / det;

  return C;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

