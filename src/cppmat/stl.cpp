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
inline std::vector<X> sort_pmodulo(const std::vector<X> &in, X n, bool reverse)
{
  std::vector<X> out = in;

  // take the modulo (e.g. to correct for 'periodicity')
  for ( auto &i : out )
    i = (n + (i%n)) % n;

  // sort
  std::sort(out.begin(),out.end());

  // reverse order
  if ( reverse ) std::reverse(out.begin(), out.end());

  return out;
}

// =================================================================================================

template<class X>
inline std::string to_string(const std::vector<X> &A)
{
  std::string out = "(";

  for ( size_t i = 0 ; i < A.size()-1 ; ++i )
    out += std::to_string(A[i]) + ", ";

  out += std::to_string(A[A.size()-1]) + ")";

  return out;
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

