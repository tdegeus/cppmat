/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_CPP
#define CPPMAT_VAR_CARTESIAN_CPP

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
cppmat::cartesian::tensor4<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor4<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          for ( size_t m = 0 ; m < ND ; ++m )
            for ( size_t n = 0 ; n < ND ; ++n )
              C(i,j,m,n) += A(i,j,k,l) * B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j) += A(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j) += A(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j) += A(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor4<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(k,l) += A(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor4<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(k,l) += A(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> ddot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor4<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(k,l) += A[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += ( A[i*ND-(i-1)*i/2] * B[i*ND-(i-1)*i/2] );

  for ( size_t i = 0 ; i<ND ; ++i )
    for ( size_t j = i+1 ; j<ND ; ++j )
      C += ( static_cast<X>(2) * A[i*ND-(i-1)*i/2+j-i] * B[i*ND-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i*ND-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B[i*ND-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X ddot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B[i];

  return C;
}

// =================================================================================================
// tensor products: dot
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += A(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += A(i,j) * B[j];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += A[i] * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C = cppmat::cartesian::tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += A[i] * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2d<X> dot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2d<X> C = cppmat::cartesian::tensor2d<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C[i] += A[i] * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i) += A(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i) += A(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += A[i] * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += A(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += A(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::vector<X> dot(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::vector<X> C = cppmat::cartesian::vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += A(i) * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X dot(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i) * B(i);

  return C;
}

// =================================================================================================
// tensor products: dyadic
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j,k,k) += A(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2s<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j,k,k) += A(i,j) * B[k];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += A[i] * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2s<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += A[i] * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> dyadic(
  const cppmat::cartesian::tensor2d<X> &A, const cppmat::cartesian::tensor2d<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C = cppmat::cartesian::tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,i,k,k) += A[i] * B[k];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> dyadic(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) = A(i) * B(j);

  return C;
}

// =================================================================================================
// cross (outer) product
// =================================================================================================

template<class X>
inline
cppmat::cartesian::vector<X> cross(
  const cppmat::cartesian::vector<X> &A, const cppmat::cartesian::vector<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  if ( ND != 3 )
    throw std::runtime_error("'cross' only implemented in 3D");

  cppmat::cartesian::vector<X> C(3);

  C[0] =                     A[1]*B[2]-B[1]*A[2] ;
  C[1] = static_cast<X>(-1)*(A[0]*B[2]-B[0]*A[2]);
  C[2] =                     A[0]*B[1]-B[0]*A[1] ;

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor4<X> T(const cppmat::cartesian::tensor4<X> &A)
{
  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(l,k,j,i) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> RT(const cppmat::cartesian::tensor4<X> &A)
{
  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,l,k) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor4<X> LT(const cppmat::cartesian::tensor4<X> &A)
{
  size_t ND = A.ndim();

  cppmat::cartesian::tensor4<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(j,i,k,l) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2<X> T(const cppmat::cartesian::tensor2<X> &A)
{
  size_t ND = A.ndim();

  cppmat::cartesian::tensor2<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j,i) = A(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2s<X> T(const cppmat::cartesian::tensor2s<X> &A)
{
  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2d<X> T(const cppmat::cartesian::tensor2d<X> &A)
{
  return A;
}

// =================================================================================================
// trace
// =================================================================================================

template<class X>
inline
X trace(const cppmat::cartesian::tensor2<X> &A)
{
  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X trace(const cppmat::cartesian::tensor2s<X> &A)
{
  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X trace(const cppmat::cartesian::tensor2d<X> &A)
{
  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i];

  return C;
}

// =================================================================================================
// determinant
// =================================================================================================

template<class X>
inline
X det(const cppmat::cartesian::tensor2<X> &A)
{
  size_t ND = A.ndim();

  if ( ND == 2 )
   return A[0] * A[3] - A[1] * A[2];

  if ( ND == 3 )
    return ( A[0] * A[4] * A[8] +
             A[1] * A[5] * A[6] +
             A[2] * A[3] * A[7] ) -
           ( A[2] * A[4] * A[6] +
             A[1] * A[3] * A[8] +
             A[0] * A[5] * A[7] );

  throw std::runtime_error("'det' only implemented in 2D/3D");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X det(const cppmat::cartesian::tensor2s<X> &A)
{
  size_t ND = A.ndim();

  if ( ND == 2 )
   return A[0] * A[2] - A[1] * A[1];

  if ( ND == 3 )
    return (                     A[0] * A[3] * A[5] +
             static_cast<X>(2) * A[1] * A[2] * A[4] ) -
           (                     A[4] * A[4] * A[0] +
                                 A[2] * A[2] * A[3] +
                                 A[1] * A[1] * A[5] );

  throw std::runtime_error("'det' only implemented in 2D/3D");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X det(const cppmat::cartesian::tensor2d<X> &A)
{
  size_t ND = A.ndim();

  X C = static_cast<X>(1);

  for ( size_t i = 0 ; i < ND ; ++i )
    C *= A[i];

  return C;
}

// =================================================================================================
// inverse
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor2<X> inv(const cppmat::cartesian::tensor2<X> &A)
{
  size_t ND = A.ndim();

  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::cartesian::tensor2<X> C(ND);

  if ( ND == 2 )
  {
    C[0] =                      A[3] / det;
    C[1] = static_cast<X>(-1) * A[1] / det;
    C[2] = static_cast<X>(-1) * A[2] / det;
    C[3] =                      A[0] / det;

    return C;
  }

  if ( ND == 3 )
  {
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

  throw std::runtime_error("'inv' only implemented in 2D/3D");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2s<X> inv(const cppmat::cartesian::tensor2s<X> &A)
{
  size_t ND = A.ndim();

  // compute determinant
  X det = A.det();

  // allocate result
  cppmat::cartesian::tensor2s<X> C(ND);

  if ( ND==2 )
  {
    C[0] =                      A[2] / det;
    C[1] = static_cast<X>(-1) * A[1] / det;
    C[2] =                      A[0] / det;

    return C;
  }

  if ( ND==3 )
  {
    C[0] = (A[3]*A[5]-A[4]*A[4]) / det;
    C[1] = (A[2]*A[4]-A[1]*A[5]) / det;
    C[2] = (A[1]*A[4]-A[2]*A[3]) / det;
    C[3] = (A[0]*A[5]-A[2]*A[2]) / det;
    C[4] = (A[2]*A[1]-A[0]*A[4]) / det;
    C[5] = (A[0]*A[3]-A[1]*A[1]) / det;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
cppmat::cartesian::tensor2d<X> inv(const cppmat::cartesian::tensor2d<X> &A)
{
  size_t ND = A.ndim();

  cppmat::cartesian::tensor2d<X> C(ND);

  for ( size_t i = 0; i < ND ; ++i )
    C[i] = static_cast<X>(1) / A[i];

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X length(const cppmat::cartesian::vector<X> &A)
{
  size_t ND = A.ndim();

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += std::pow(A[i],2.);

  return std::sqrt(C);
}

// =================================================================================================

}} // namespace ...

#endif

