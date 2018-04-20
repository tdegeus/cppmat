/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_STL_CPP
#define CPPMAT_STL_CPP

// -------------------------------------------------------------------------------------------------

#include "stl.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================

template<class X>
inline std::vector<X> del(const std::vector<X> &A, int idx)
{
  int n = static_cast<int>(A.size());

  idx = ( idx < 0 ) ? idx + n : ( idx >= n ) ? idx - n : idx ;

  assert( idx >= 0 );
  assert( idx  < n );

  std::vector<X> B = A;

  B.erase(B.begin()+idx, B.begin()+idx+1);

  return B;
}

// =================================================================================================

template<class X>
inline std::vector<X> del(const std::vector<X> &A, size_t idx)
{
  assert( idx < A.size() );

  std::vector<X> B = A;

  B.erase(B.begin()+idx, B.begin()+idx+1);

  return B;
}

// =================================================================================================

template<class X>
inline X abs(X A)
{
  return std::abs(A);
}

// =================================================================================================

template<class X>
inline matrix<X> abs(const matrix<X> &A)
{
  matrix<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X>
inline cartesian2d::tensor2<X> abs(const cartesian2d::tensor2<X> &A)
{
  cartesian2d::tensor2<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X>
inline cartesian2d::tensor2s<X> abs(const cartesian2d::tensor2s<X> &A)
{
  cartesian2d::tensor2s<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

} // namespace ...

#endif

