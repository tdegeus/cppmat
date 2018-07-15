/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PRIVATE_HPP
#define CPPMAT_PRIVATE_HPP

// -------------------------------------------------------------------------------------------------

#include "private.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace Private {

// =================================================================================================

template<class X>
inline
std::vector<X> sort_axes(const std::vector<X> &in, X n, bool reverse)
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

inline bool equal(double a, double b)
{
  return std::fabs(a - b) <= std::numeric_limits<double>::epsilon();
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

