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
// tensor products
// =================================================================================================

template<class X>
inline
cppmat::cartesian::tensor4<X> ddot(
  const cppmat::cartesian::tensor4<X> &A, const cppmat::cartesian::tensor4<X> &B
)
{
  assert( A.ndim() == B.ndim() );

  size_t ND = A.ndim();

  tensor4<X> C = tensor4<X>::Zero(ND);

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

  tensor2<X> C = tensor2<X>::Zero(ND);

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

  tensor2<X> C = tensor2<X>::Zero(ND);

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

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j) += A(i,j,k,k) * B(k,k);

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

  tensor4<X> C(ND);

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

  tensor4<X> C(ND);

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

  tensor4<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(j,i,k,l) = A(i,j,k,l);

  return C;
}

// =================================================================================================

}} // namespace ...

#endif

