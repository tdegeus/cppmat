/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VECTOR_CPP
#define CPPMAT_VECTOR_CPP

// -------------------------------------------------------------------------------------------------

#include "regular_vector.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline vector<X>::vector(size_t n)
{
  // store shape, and other size parameters, allocate "mData"
  resize(n);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Arange(size_t n)
{
  // call basic constructor
  vector<X> out(n);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Zero(size_t n)
{
  // call basic constructor
  vector<X> out(n);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Ones(size_t n)
{
  // call basic constructor
  vector<X> out(n);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Constant(size_t n, X D)
{
  // call basic constructor
  vector<X> out(n);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline vector<X> vector<X>::Copy(Iterator first, Iterator last)
{
  // call basic constructor
  vector<X> out(last-first);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void vector<X>::resize(size_t n)
{
  // update size
  N     = n;
  mSize = n;

  // allocate data
  mData.resize(mSize);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t vector<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::ndim() const
{
  return 1;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::shape(int i) const
{
  // check axis: (0) or (-1)
  assert( i  <  1 );
  assert( i >= -1 );

  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::shape(size_t i) const
{
  // check axis: (0)
  assert( i < 1 );

  return N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> ret(1);

  ret[0] = N;

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::strides(bool bytes) const
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

template<class X>
inline X& vector<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline X& vector<X>::operator()(size_t a)
{
  assert( a < N );

  return mData[a];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator()(size_t a) const
{
  assert( a < N );

  return mData[a];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline X* vector<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* vector<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline auto vector<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto vector<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline auto vector<X>::item(size_t a)
{
  assert( a < N );

  return begin() + a;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::item(size_t a) const
{
  assert( a < N );

  return begin() + a;
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void vector<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void vector<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last - first );

  std::copy(first, last, begin());
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline vector<X>& vector<X>::operator*= (const vector<X> &B)
{
  assert( N == B.shape(0) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  assert( N == B.shape(0) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  assert( N == B.shape(0) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const vector<X> &B)
{
  assert( N == B.shape(0) );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  assert( A.shape(0) == B.shape(0) );

  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C(A.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C(B.shape(0));

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X>
inline size_t vector<X>::argmin() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::argmax() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X>
inline X vector<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X>
inline X vector<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X>
inline X vector<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X>
inline double vector<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X>
inline double vector<X>::average(const vector<X> &weights, bool norm) const
{
  assert( N == weights.shape(0) );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// =================================================================================================
// find all non-zero entries
// =================================================================================================

template<class X>
inline vector<size_t> vector<X>::where() const
{
  size_t nnz = 0;

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] )
      ++nnz;

  vector<size_t> out(nnz);

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
// basic algebra : absolute value
// =================================================================================================

template<class X>
inline void vector<X>::abs()
{
  for ( auto &i : mData )
    i = std::abs(i);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < N ; ++j ) {
    if ( j != N-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else            std::printf((fmt + ";\n").c_str(), (*this)(j));
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < src.size() ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src(j);
    if ( j != src.size()-1 ) out << ", ";
  }

  return out;
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

