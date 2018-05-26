/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_3_CPP
#define CPPMAT_FIX_CARTESIAN_3_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// tensor products: ddot
// =================================================================================================

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2<X,3> &A, const cppmat::tiny::cartesian::tensor2d<X,3> &B
)
{
  return A[0]*B[0] + A[4]*B[1] + A[8]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,3> &A, const cppmat::tiny::cartesian::tensor2s<X,3> &B
)
{
  X C;

  C  = A[0] * B[0];
  C += A[1] * B[1] * static_cast<X>(2);
  C += A[2] * B[2] * static_cast<X>(2);
  C += A[3] * B[3];
  C += A[4] * B[4] * static_cast<X>(2);
  C += A[5] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,3> &A, const cppmat::tiny::cartesian::tensor2d<X,3> &B
)
{
  return A[0]*B[0] + A[3]*B[1] + A[5]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,3> &A, const cppmat::tiny::cartesian::tensor2<X,3> &B
)
{
  return A[0]*B[0] + A[1]*B[4] + A[2]*B[8];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,3> &A, const cppmat::tiny::cartesian::tensor2s<X,3> &B
)
{
  return A[0]*B[0] + A[1]*B[3] + A[2]*B[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,3> &A, const cppmat::tiny::cartesian::tensor2d<X,3> &B
)
{
  return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

// =================================================================================================
// cross (outer) product
// =================================================================================================

template<class X>
inline
cppmat::tiny::cartesian::vector<X,3> cross(
  const cppmat::tiny::cartesian::vector<X,3> &A, const cppmat::tiny::cartesian::vector<X,3> &B
)
{
  cppmat::tiny::cartesian::vector<X,3> C;

  C[0] =                     A[1]*B[2]-B[1]*A[2] ;
  C[1] = static_cast<X>(-1)*(A[0]*B[2]-B[0]*A[2]);
  C[2] =                     A[0]*B[1]-B[0]*A[1] ;

  return C;
}

// =================================================================================================
// determinant
// =================================================================================================

template<class X>
inline
X det(const cppmat::tiny::cartesian::tensor2<X,3> &A)
{
  return ( A[0] * A[4] * A[8] +
           A[1] * A[5] * A[6] +
           A[2] * A[3] * A[7] ) -
         ( A[2] * A[4] * A[6] +
           A[1] * A[3] * A[8] +
           A[0] * A[5] * A[7] );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X det(const cppmat::tiny::cartesian::tensor2s<X,3> &A)
{
  return (                     A[0] * A[3] * A[5] +
           static_cast<X>(2) * A[1] * A[2] * A[4] ) -
         (                     A[4] * A[4] * A[0] +
                               A[2] * A[2] * A[3] +
                               A[1] * A[1] * A[5] );
}

// =================================================================================================
// inverse
// =================================================================================================

template<class X>
inline
cppmat::tiny::cartesian::tensor2<X,3> inv(const cppmat::tiny::cartesian::tensor2<X,3> &A)
{
  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::tiny::cartesian::tensor2<X,3> C;

  // compute inverse
  C[0] = (A[4]*A[8]-A[5]*A[7]) / det;
  C[1] = (A[2]*A[7]-A[1]*A[8]) / det;
  C[2] = (A[1]*A[5]-A[2]*A[4]) / det;
  C[3] = (A[5]*A[6]-A[3]*A[8]) / det;
  C[4] = (A[0]*A[8]-A[2]*A[6]) / det;
  C[5] = (A[2]*A[3]-A[0]*A[5]) / det;
  C[6] = (A[3]*A[7]-A[4]*A[6]) / det;
  C[7] = (A[1]*A[6]-A[0]*A[7]) / det;
  C[8] = (A[0]*A[4]-A[1]*A[3]) / det;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::tiny::cartesian::tensor2s<X,3> inv(const cppmat::tiny::cartesian::tensor2s<X,3> &A)
{
  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::tiny::cartesian::tensor2s<X,3> C;

  // compute inverse
  C[0] = (A[3]*A[5]-A[4]*A[4]) / det;
  C[1] = (A[2]*A[4]-A[1]*A[5]) / det;
  C[2] = (A[1]*A[4]-A[2]*A[3]) / det;
  C[3] = (A[0]*A[5]-A[2]*A[2]) / det;
  C[4] = (A[2]*A[1]-A[0]*A[4]) / det;
  C[5] = (A[0]*A[3]-A[1]*A[1]) / det;

  return C;
}

// =================================================================================================

}} // namespace ...

#endif

