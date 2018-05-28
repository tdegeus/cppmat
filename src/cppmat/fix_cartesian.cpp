/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_CPP
#define CPPMAT_FIX_CARTESIAN_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// tensor products: ddot
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor4<X,ND> &A, const cppmat::tiny::cartesian::tensor4<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

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

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor4<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j) += A(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor4<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j) += A(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor4<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j) += A(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor4<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(k,l) += A(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor4<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(k,l) += A(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> ddot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor4<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(k,l) += A[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += A(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += ( A[i*ND-(i-1)*i/2] * B[i*ND-(i-1)*i/2] );

  for ( size_t i = 0 ; i<ND ; ++i )
    for ( size_t j = i+1 ; j<ND ; ++j )
      C += ( static_cast<X>(2) * A[i*ND-(i-1)*i/2+j-i] * B[i*ND-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i*ND-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B[i*ND-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X ddot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A[i]*B[i];

  return C;
}

// =================================================================================================
// tensor products: dot
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += A(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += A(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += A(i,j) * B[j];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += A[i] * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C = cppmat::tiny::cartesian::tensor2<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += A[i] * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2d<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2d<X,ND> C = cppmat::tiny::cartesian::tensor2d<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    C[i] += A[i] * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i) += A(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i) += A(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += A[i] * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += A(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += A(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> dot(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::vector<X,ND> C = cppmat::tiny::cartesian::vector<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += A(i) * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X dot(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += A(i) * B(i);

  return C;
}

// =================================================================================================
// tensor products: dyadic
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j,k,k) += A(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += A(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2s<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j,k,k) += A(i,j) * B[k];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += A[i] * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2s<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += A[i] * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> dyadic(
  const cppmat::tiny::cartesian::tensor2d<X,ND> &A, const cppmat::tiny::cartesian::tensor2d<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C = cppmat::tiny::cartesian::tensor4<X,ND>::Zero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,i,k,k) += A[i] * B[k];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> dyadic(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C;

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) = A(i) * B(j);

  return C;
}

// =================================================================================================
// cross (outer) product
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::vector<X,ND> cross(
  const cppmat::tiny::cartesian::vector<X,ND> &A, const cppmat::tiny::cartesian::vector<X,ND> &B
)
{
  assert( false );

  UNUSED(B);

  return A;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> T(const cppmat::tiny::cartesian::tensor4<X,ND> &A)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C;

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(l,k,j,i) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> RT(const cppmat::tiny::cartesian::tensor4<X,ND> &A)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C;

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,l,k) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor4<X,ND> LT(const cppmat::tiny::cartesian::tensor4<X,ND> &A)
{
  cppmat::tiny::cartesian::tensor4<X,ND> C;

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(j,i,k,l) = A(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> T(const cppmat::tiny::cartesian::tensor2<X,ND> &A)
{
  cppmat::tiny::cartesian::tensor2<X,ND> C;

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j,i) = A(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2s<X,ND> T(const cppmat::tiny::cartesian::tensor2s<X,ND> &A)
{
  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2d<X,ND> T(const cppmat::tiny::cartesian::tensor2d<X,ND> &A)
{
  return A;
}

// =================================================================================================
// trace
// =================================================================================================

template<class X, size_t ND>
inline
X trace(const cppmat::tiny::cartesian::tensor2<X,ND> &A)
{
  X C = A[0];

  for ( size_t i = 1 ; i < ND ; ++i )
    C += A(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X trace(const cppmat::tiny::cartesian::tensor2s<X,ND> &A)
{
  X C = A[0];

  for ( size_t i = 1 ; i < ND ; ++i )
    C += A(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X trace(const cppmat::tiny::cartesian::tensor2d<X,ND> &A)
{
  X C = A[0];

  for ( size_t i = 1 ; i < ND ; ++i )
    C += A[i];

  return C;
}

// =================================================================================================
// determinant
// =================================================================================================

template<class X, size_t ND>
inline
X det(const cppmat::tiny::cartesian::tensor2<X,ND> &A)
{
  assert( false );

  UNUSED(A);

  return static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X det(const cppmat::tiny::cartesian::tensor2s<X,ND> &A)
{
  assert( false );

  UNUSED(A);

  return static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X det(const cppmat::tiny::cartesian::tensor2d<X,ND> &A)
{
  X C = A[0];

  for ( size_t i = 1 ; i < ND ; ++i )
    C *= A[i];

  return C;
}

// =================================================================================================
// inverse
// =================================================================================================

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2<X,ND> inv(const cppmat::tiny::cartesian::tensor2<X,ND> &A)
{
  assert(false);

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2s<X,ND> inv(const cppmat::tiny::cartesian::tensor2s<X,ND> &A)
{
  assert(false);

  return A;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
cppmat::tiny::cartesian::tensor2d<X,ND> inv(const cppmat::tiny::cartesian::tensor2d<X,ND> &A)
{
  cppmat::tiny::cartesian::tensor2d<X,ND> C;

  for ( size_t i = 0; i < ND ; ++i )
    C[i] = static_cast<X>(1) / A[i];

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X, size_t ND>
inline
X length(const cppmat::tiny::cartesian::vector<X,ND> &A)
{
  X C = std::pow(A[0],2.);

  for ( size_t i = 1 ; i < ND ; ++i )
    C += std::pow(A[i],2.);

  return std::sqrt(C);
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

