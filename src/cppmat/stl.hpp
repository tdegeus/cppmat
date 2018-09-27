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

template<typename X>
inline
std::vector<X> del(const std::vector<X> &A, int idx)
{
  int n = static_cast<int>(A.size());

  idx = ( idx < 0 ) ? idx + n : ( idx >= n ) ? idx - n : idx ;

  Assert( idx >= 0 );
  Assert( idx  < n );

  std::vector<X> B = A;

  B.erase(B.begin()+idx, B.begin()+idx+1);

  return B;
}

// =================================================================================================
// delete a specific item from a vector
// =================================================================================================

template<typename X>
inline
std::vector<X> del(const std::vector<X> &A, size_t idx)
{
  Assert( idx < A.size() );

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

template<typename X>
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
// linearly spaced array
// =================================================================================================

template <typename T>
std::vector<T> linspace(T a, T b, size_t N)
{
  // spacing
  double h = (static_cast<double>(b) - static_cast<double>(a)) / static_cast<double>(N-1);

  // allocate output
  std::vector<double> out(N);

  // temporary variables
  typename std::vector<double>::iterator x;
  double val;

  // loop to fill
  for ( x = out.begin(), val = a; x != out.end(); ++x, val += h)
    *x = val;

  std::vector<T> outT(out.begin(), out.end());

  return outT;
}

// =================================================================================================
// minimum/maximum from a vector
// =================================================================================================

template<typename X> X min(const std::vector<X> &A)
{
  return *std::min_element(A.begin(),A.end());
}

// -------------------------------------------------------------------------------------------------

template<typename X> X max(const std::vector<X> &A)
{
  return *std::max_element(A.begin(),A.end());
}

// =================================================================================================
// minimum/maximum from two vectors of equal size
// =================================================================================================

template<typename X> std::vector<X> min(const std::vector<X> &A, const std::vector<X> &B)
{
  Assert( A.size () == B.size() );

  std::vector<X> C(A.size());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = std::min(A[i], B[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<typename X> std::vector<X> max(const std::vector<X> &A, const std::vector<X> &B)
{
  Assert( A.size () == B.size() );

  std::vector<X> C(A.size());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = std::max(A[i], B[i]);

  return C;
}

// =================================================================================================

} // namespace ...

// =================================================================================================
// print operator
// =================================================================================================

#ifndef CPPMAT_NOSTD
template<typename X>
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

