/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_STL_HPP
#define CPPMAT_STL_HPP

// -------------------------------------------------------------------------------------------------

#include "stl.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// delete a specific item from a vector
// =================================================================================================

template<class X>
inline
std::vector<X> del(const std::vector<X> &A, int idx)
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
// delete a specific item from a vector
// =================================================================================================

template<class X>
inline
std::vector<X> del(const std::vector<X> &A, size_t idx)
{
  assert( idx < A.size() );

  std::vector<X> B = A;

  B.erase(B.begin()+idx, B.begin()+idx+1);

  return B;
}

// =================================================================================================
// return the indices that would sort the vector
// =================================================================================================

template <typename X>
inline
std::vector<size_t> argsort(const std::vector<X> &v, bool ascending)
{
  // initialize original index locations
  // - allocate
  std::vector<size_t> idx(v.size());
  // - fill
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in "v"
  std::sort(idx.begin(), idx.end(), [&v](size_t i, size_t j) {return v[i] < v[j];});

  // no reversing of order
  if ( ascending ) return idx;

  // reverse order
  // - allocate
  std::vector<size_t> jdx(v.size());
  // - fill
  for ( size_t i = 0 ; i < v.size() ; ++i ) jdx[v.size()-i-1] = idx[i];
  // - return
  return jdx;
}

// =================================================================================================
// join items to string
// =================================================================================================

template<class X>
inline
std::string to_string(const std::vector<X> &A, std::string join)
{
  std::string out = "";

  for ( size_t i = 0 ; i < A.size()-1 ; ++i )
    out += std::to_string(A[i]) + join;

  out += std::to_string(A[A.size()-1]);

  return out;
}

// =================================================================================================

} // namespace ...

// =================================================================================================
// print operator
// =================================================================================================

#ifndef CPPMAT_NOSTD
template<class X>
inline
std::ostream& operator<<(std::ostream& out, const std::vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < src.size() ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src[j];
    if ( j != src.size()-1 ) out << ", ";
  }

  return out;
}
#endif

// =================================================================================================

#endif

