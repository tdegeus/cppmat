/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_ARRAY_CPP
#define CPPMAT_ARRAY_CPP

// -------------------------------------------------------------------------------------------------

#include "regular_array.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline array<X>::array(const std::vector<size_t> &shape)
{
  // store shape, and other size parameters, allocate "mData"
  resize(shape);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::Arange(const std::vector<size_t> &shape)
{
  // call basic constructor
  array<X> out(shape);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::Zero(const std::vector<size_t> &shape)
{
  // call basic constructor
  array<X> out(shape);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::Ones(const std::vector<size_t> &shape)
{
  // call basic constructor
  array<X> out(shape);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::Constant(const std::vector<size_t> &shape, X D)
{
  // call basic constructor
  array<X> out(shape);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline array<X> array<X>::Copy(const std::vector<size_t> &shape, Iterator first, Iterator last)
{
  // call basic constructor
  array<X> out(shape);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void array<X>::resize(const std::vector<size_t> &shape)
{
  assert( shape.size()  > 0       );
  assert( shape.size() <= MAX_DIM );

  // update number of dimensions
  mNdim = shape.size();

  // initialize shape in all directions to "1"
  for ( size_t i = 0 ; i < MAX_DIM ; ++i )
  {
    mShape  [i] = 1;
    mStrides[i] = 1;
  }

  // initialize size
  mSize = 1;

  // update shape/size
  for ( size_t i = 0 ; i < mNdim ; ++i )
  {
    mShape[i] = shape[i];
    mSize    *= shape[i];
  }

  // set storage strides
  for ( size_t i = 0 ; i < mNdim ; ++i )
    for ( size_t j = i+1 ; j < mNdim ; ++j )
      mStrides[i] *= mShape[j];

  // allocate data
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::reshape(const std::vector<size_t> &shape)
{
  // check that the size is unchanged
  #ifndef NDEBUG
    size_t n = 1;

    for ( auto &i : shape ) n *= i;

    assert( n == mSize );
  #endif

  // process new shape
  resize(shape);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::chdim(size_t ndim)
{
  // check that all removed dimensions are of shape 1
  #ifndef NDEBUG
    for ( size_t i = ndim ; i < MAX_DIM ; ++i ) assert( mShape[i] == 1 );
  #endif

  // update number of dimensions
  mNdim = ndim;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t array<X>::size() const
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::ndim() const
{
  return mNdim;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::shape(int i) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( i  <      static_cast<int>(mNdim) );
  assert( i >= -1 * static_cast<int>(mNdim) );

  // get number of dimensions as integer
  int n = static_cast<int>(mNdim);

  // correct periodic index
  i = ( n + (i%n) ) % n;

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::shape(size_t i) const
{
  // check axis: (0,1,...,ndim-1)
  assert( i < mNdim );

  // return shape
  return mShape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> array<X>::shape() const
{
  std::vector<size_t> ret(mNdim);

  std::copy(std::begin(mShape), std::begin(mShape)+mNdim, ret.begin());

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> array<X>::strides(bool bytes) const
{
  std::vector<size_t> ret(mNdim);

  std::copy(std::begin(mStrides), std::begin(mStrides)+mNdim, ret.begin());

  if ( bytes )
    for ( size_t i = 0 ; i < mNdim ; ++i )
      ret[i] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X>
inline X& array<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline X& array<X>::operator()(size_t a)
{
  assert( mNdim >= 1 );

  assert( a < mShape[0] );

  return mData[a*mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a) const
{
  assert( mNdim >= 1 );

  assert( a < mShape[0] );

  return mData[a*mStrides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b)
{
  assert( mNdim >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[a*mStrides[0]+b*mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b) const
{
  assert( mNdim >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return mData[a*mStrides[0]+b*mStrides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c)
{
  assert( mNdim >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c) const
{
  assert( mNdim >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d)
{
  assert( mNdim >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mNdim >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( mNdim >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mNdim >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( mNdim >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return \
  mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]+f*mStrides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(
  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mNdim >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return \
  mData[a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]+f*mStrides[5]];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X>
template<class Iterator>
inline X& array<X>::at(Iterator first, Iterator last)
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0     );
  assert( static_cast<size_t>(last-first) <= mNdim );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    assert( (*it) < (*shape) );
    // - update index
    idx += (*it) * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return mData[idx];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline const X& array<X>::at(Iterator first, Iterator last) const
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0     );
  assert( static_cast<size_t>(last-first) <= mNdim );

  // iterator to shape and stride
  size_t *shape  = &mShape  [0];
  size_t *stride = &mStrides[0];

  // zero-initialize plain storage index
  size_t idx = 0;

  // loop over array-indices
  for ( auto it = first ; it != last ; ++it )
  {
    // - check array index
    assert( (*it) < (*shape) );
    // - update index
    idx += (*it) * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return mData[idx];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline size_t array<X>::compress(size_t a) const
{
  assert( mNdim >= 1 );

  assert( a < mShape[0] );

  return a*mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b) const
{
  assert( mNdim >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return a*mStrides[0]+b*mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c) const
{
  assert( mNdim >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return a*mStrides[0]+b*mStrides[1]+c*mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mNdim >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mNdim >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mNdim >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]+f*mStrides[5];
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X>
inline std::vector<size_t> array<X>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(mNdim);

  // reconstruct
  for ( size_t j = 0 ; j < mNdim ; ++j ) {
    idx[j] = (i - i%mStrides[j]) / mStrides[j];
    i -= idx[j] * mStrides[j];
  }

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline X* array<X>::data()
{
  return mData.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* array<X>::data() const
{
  return mData.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline auto array<X>::begin()
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::begin() const
{
  return mData.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::end()
{
  return mData.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::end() const
{
  return mData.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto array<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline auto array<X>::item(size_t a)
{
  assert( mNdim >= 1 );

  assert( a < mShape[0] );

  return begin() + a*mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a) const
{
  assert( mNdim >= 1 );

  assert( a < mShape[0] );

  return begin() + a*mStrides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b)
{
  assert( mNdim >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() + a*mStrides[0]+b*mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b) const
{
  assert( mNdim >= 2 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );

  return begin() + a*mStrides[0]+b*mStrides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c)
{
  assert( mNdim >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return begin() + a*mStrides[0]+b*mStrides[1]+c*mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c) const
{
  assert( mNdim >= 3 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );

  return begin() + a*mStrides[0]+b*mStrides[1]+c*mStrides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d)
{
  assert( mNdim >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return begin() + a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d) const
{
  assert( mNdim >= 4 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );

  return begin() + a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( mNdim >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return begin()+a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( mNdim >= 5 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );

  return begin()+a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( mNdim >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return begin()+a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]+f*mStrides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( mNdim >= 6 );

  assert( a < mShape[0] );
  assert( b < mShape[1] );
  assert( c < mShape[2] );
  assert( d < mShape[3] );
  assert( e < mShape[4] );
  assert( f < mShape[5] );

  return begin()+a*mStrides[0]+b*mStrides[1]+c*mStrides[2]+d*mStrides[3]+e*mStrides[4]+f*mStrides[5];
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void array<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void array<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last - first );

  std::copy(first, last, begin());
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline array<X>& array<X>::operator*= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator/= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator+= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator-= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator* (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.ndim()  == B.ndim()  );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator/ (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.ndim()  == B.ndim()  );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator+ (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.ndim()  == B.ndim()  );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator- (const array<X> &A, const array<X> &B)
{
  assert( A.shape() == B.shape() );
  assert( A.ndim()  == B.ndim()  );

  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator* (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator/ (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator+ (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator- (const array<X> &A, const X &B)
{
  array<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator* (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator/ (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator+ (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> operator- (const X &A, const array<X> &B)
{
  array<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra : location of the minimum/maximum
// =================================================================================================

template<class X>
inline std::vector<size_t> array<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> array<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra : minimum
// =================================================================================================

template<class X>
inline X array<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::minCoeff(size_t axis) const
{
  // check input
  assert( axis < mNdim );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(this->shape(),axis), this->maxCoeff());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = this->strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::min(out[ni], mData[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::minCoeff(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(mNdim) );
  assert( axis >= -1 * static_cast<int>(mNdim) );

  // get number of dimensions as integer
  int n = static_cast<int>(mNdim);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return this->minCoeff(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::minCoeff(const std::vector<int> &axes_in) const
{
  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(mNdim), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.minCoeff(axis);

  return out;
}

// =================================================================================================
// basic algebra : maximum
// =================================================================================================

template<class X>
inline X array<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::maxCoeff(size_t axis) const
{
  // check input
  assert( axis < mNdim );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(this->shape(),axis), this->minCoeff());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = this->strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::max(out[ni], mData[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::maxCoeff(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(mNdim) );
  assert( axis >= -1 * static_cast<int>(mNdim) );

  // get number of dimensions as integer
  int n = static_cast<int>(mNdim);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return this->maxCoeff(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::maxCoeff(const std::vector<int> &axes_in) const
{
  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(mNdim), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.maxCoeff(axis);

  return out;
}

// =================================================================================================
// basic algebra : sum
// =================================================================================================

template<class X>
inline X array<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::sum(size_t axis) const
{
  // zero-initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Zero(del(this->shape(),axis));

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = this->strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), mSize);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < mSize ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] += mData[i];
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::sum(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(mNdim) );
  assert( axis >= -1 * static_cast<int>(mNdim) );

  // get number of dimensions as integer
  int n = static_cast<int>(mNdim);

  // correct periodic axis
  axis = ( n + (axis%n) ) % n;

  // compute
  return this->sum(static_cast<size_t>(axis));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::sum(const std::vector<int> &axes_in) const
{
  // correct for 'periodicity', sort from high to low
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(mNdim), true);

  // copy array
  array<X> out = (*this);

  // loop to compute
  for ( auto &axis : axes )
    out = out.sum(axis);

  return out;
}

// =================================================================================================
// basic algebra : mean
// =================================================================================================

template<class X>
inline double array<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::mean(size_t axis) const
{
  array<X> weights = array<X>::Ones(this->shape());

  return (*this).sum(axis) / weights.sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::mean(int axis) const
{
  array<X> weights = array<X>::Ones(this->shape());

  return (*this).sum(axis) / weights.sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::mean(const std::vector<int> &axes) const
{
  array<X> weights = array<X>::Ones(this->shape());

  return (*this).sum(axes) / weights.sum(axes);
}

// =================================================================================================
// basic algebra : weighted average
// =================================================================================================

template<class X>
inline double array<X>::average(const array<X> &weights, bool norm) const
{
  assert( shape() == weights.shape() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::average(const array<X> &weights, size_t axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::average(const array<X> &weights, int axis, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axis) / weights.sum(axis);
  else        return (weights*(*this)).sum(axis);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::average(
  const array<X> &weights, const std::vector<int> &axes, bool norm) const
{
  if ( norm ) return (weights*(*this)).sum(axes) / weights.sum(axes);
  else        return (weights*(*this)).sum(axes);
}

// =================================================================================================
// basic algebra : absolute value
// =================================================================================================

template<class X>
inline array<X> array<X>::abs() const
{
  array<X> out = (*this);

  for ( auto &i : out )
    i = std::abs(i);

  return out;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void array<X>::printf(std::string fmt) const
{
  if ( mNdim == 1 )
  {
    for ( size_t j = 0 ; j < mShape[0] ; ++j ) {
      if ( j != mShape[0]-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
      else                    std::printf((fmt + ";\n").c_str(), (*this)(j));
    }

    return;
  }

  if ( mNdim == 2 )
  {
    for ( size_t i = 0 ; i < mShape[0] ; ++i ) {
      for ( size_t j = 0 ; j < mShape[1] ; ++j ) {
        if ( j != mShape[1]-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
        else                    std::printf((fmt + ";\n").c_str(), (*this)(i,j));
      }
    }

    return;
  }

  std::cout << "cppmat::array[";

  for ( size_t i = 0 ; i < mNdim-1 ; ++i )
    std::cout << shape(i) << ",";

  std::cout << shape(mNdim-1) << "]\n";

}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const array<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  if ( src.ndim() == 1 )
  {
    for ( size_t j = 0 ; j < src.shape(0) ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(j);
      if ( j != src.shape(0)-1 ) out << ", ";
    }

    return out;
  }

  if ( src.ndim() == 2 )
  {
    for ( size_t i = 0 ; i < src.shape(0) ; ++i ) {
      for ( size_t j = 0 ; j < src.shape(1) ; ++j ) {
        out << std::setw(w) << std::setprecision(p) << src(i,j);
        if      ( j != src.shape(1)-1 ) out << ", ";
        else if ( i != src.shape(0)-1 ) out << ";" << std::endl;
        else                            out << ";";
      }
    }

    return out;
  }

  out << "cppmat::array[";

  for ( size_t i = 0 ; i < src.ndim()-1 ; ++i )
    out << src.shape(i) << ",";

  out << src.shape(src.ndim()-1) << "]";

  return out;
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

