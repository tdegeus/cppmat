/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TINY_VECTOR_CPP
#define CPPMAT_VIEW_TINY_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "view_tiny_vector.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace tiny {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t N>
inline vector<X,N> vector<X,N>::Map(const X *D)
{
  // call basic constructor
  vector<X,N> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::Size()
{
  return N;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X, size_t N>
inline void vector<X,N>::setMap(const X *D)
{
  mData = D;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X, size_t N>
inline size_t vector<X,N>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::ndim() const
{
  return 1;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::shape(int i) const
{
  // check axis: (0) or (-1)
  assert( i  <  1 );
  assert( i >= -1 );

  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::shape(size_t i) const
{
  // check axis: (0)
  assert( i < 1 );

  return N;

}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline std::vector<size_t> vector<X,N>::shape() const
{
  std::vector<size_t> ret(1);

  ret[0] = N;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline std::vector<size_t> vector<X,N>::strides(bool bytes) const
{
  std::vector<size_t> ret(1);

  ret[0] = 1;

  if ( bytes )
    ret[0] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X, size_t N>
inline const X& vector<X,N>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline const X& vector<X,N>::operator()(size_t a) const
{
  assert( a < N );

  return mData[a];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X, size_t N>
inline const X* vector<X,N>::data() const
{
  return std::begin(mData);
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X, size_t N>
inline auto vector<X,N>::begin() const
{
  return std::begin(mData);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline auto vector<X,N>::end() const
{
  return std::begin(mData) + N;
}

// =================================================================================================
// iterators : index()
// =================================================================================================
template<class X, size_t N>
inline auto vector<X,N>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X, size_t N>
inline auto vector<X,N>::item(size_t a) const
{
  assert( a < N );

  return begin() + a;
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X, size_t N>
inline reg::vector<X,N> operator* (const vector<X,N> &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const vector<X,N> &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const vector<X,N> &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator- (const vector<X,N> &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator* (const vector<X,N> &A, const X &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const vector<X,N> &A, const X &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const vector<X,N> &A, const X &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator- (const vector<X,N> &A, const X &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator* (const X &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator/ (const X &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator+ (const X &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline reg::vector<X,N> operator- (const X &A, const vector<X,N> &B)
{
  reg::vector<X,N> C;

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X, size_t N>
inline size_t vector<X,N>::argmin() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::argmax() const
{
  return std::max_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline size_t vector<X,N>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X, size_t N>
inline X vector<X,N>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X, size_t N>
inline X vector<X,N>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X, size_t N>
inline X vector<X,N>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X, size_t N>
inline double vector<X,N>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X, size_t N>
inline double vector<X,N>::average(const vector<X,N> &weights, bool norm) const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// =================================================================================================
// find all non-zero entries
// =================================================================================================

template<class X, size_t N>
inline cppmat::vector<size_t> vector<X,N>::where() const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  cppmat::vector<size_t> out(nnz);

  size_t j = 0;

  for ( size_t i = 0 ; i < mSize ; ++i ) {
    if ( mData[i] ) {
      out[j] = i;
      ++j;
    }
  }

  return out;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X, size_t N>
inline void vector<X,N>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < N ; ++j ) {
    if ( j != N-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else            std::printf((fmt + ";\n").c_str(), (*this)(j));
  }
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t N>
inline std::ostream& operator<<(std::ostream& out, const vector<X,N>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < N ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src(j);
    if ( j != N-1 ) out << ", ";
  }

  return out;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

