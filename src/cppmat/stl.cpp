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
inline std::string to_string(const std::vector<X> &A)
{
  std::string out = "(";

  for ( size_t i = 0 ; i < A.size()-1 ; ++i )
    out += std::to_string(A[i]) + ", ";

  out += std::to_string(A[A.size()-1]) + ")";

  return out;
}

// =================================================================================================

template <typename X>
std::vector<size_t> argsort(const std::vector<X> &v, bool ascending) {

  // initialize original index locations
  // - allocate
  std::vector<size_t> idx(v.size());
  // - fill
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in "v"
  sort(idx.begin(), idx.end(), [&v](size_t i, size_t j) {return v[i] < v[j];});

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

template<class X>
inline array<X> abs(const array<X> &A)
{
  array<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
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
inline vector<X> abs(const vector<X> &A)
{
  vector<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X>
inline periodic::array<X> abs(const periodic::array<X> &A)
{
  periodic::array<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X>
inline periodic::matrix<X> abs(const periodic::matrix<X> &A)
{
  periodic::matrix<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X>
inline periodic::vector<X> abs(const periodic::vector<X> &A)
{
  periodic::vector<X> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X, size_t m, size_t n>
inline tiny::matrix<X,m,n> abs(const tiny::matrix<X,m,n> &A)
{
  tiny::matrix<X,m,n> out = A;

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================

template<class X, size_t n>
inline tiny::vector<X,n> abs(const tiny::vector<X,n> &A)
{
  tiny::vector<X,n> out = A;

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

