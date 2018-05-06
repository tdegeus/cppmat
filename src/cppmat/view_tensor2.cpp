/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TENSOR2_CPP
#define CPPMAT_VIEW_TENSOR2_CPP

// -------------------------------------------------------------------------------------------------

#include "view_tensor2.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian2d {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::Map(const X *D)
{
  // call basic constructor
  tensor4<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor4<X>::Size()
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Map(const X *D)
{
  // call basic constructor
  tensor2<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2<X>::Size()
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Map(const X *D)
{
  // call basic constructor
  tensor2s<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2s<X>::Size()
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d()
{
  // set dummy parameter
  mZero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Map(const X *D)
{
  // call basic constructor
  tensor2d<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2d<X>::Size()
{
  return mSize;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Map(const X *D)
{
  // call basic constructor
  vector<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::Size()
{
  return mSize;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X> inline void tensor4 <X>::setMap(const X *D) { mData = D; }
template<class X> inline void tensor2 <X>::setMap(const X *D) { mData = D; }
template<class X> inline void tensor2s<X>::setMap(const X *D) { mData = D; }
template<class X> inline void tensor2d<X>::setMap(const X *D) { mData = D; }
template<class X> inline void vector  <X>::setMap(const X *D) { mData = D; }

// =================================================================================================
// cast to different class
// =================================================================================================

template<>
template<>
inline reg::tensor2s<double> tensor2<double>::cast<reg::tensor2s<double>>() const
{
  reg::tensor2s<double> out;

  out[0] =   mData[0];
  out[1] = ( mData[1] + mData[2] ) / 2.;
  out[2] =   mData[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2d<double> tensor2<double>::cast<reg::tensor2d<double>>() const
{
  reg::tensor2d<double> out;

  out[0] = mData[0];
  out[1] = mData[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2<double> tensor2s<double>::cast<reg::tensor2<double>>() const
{
  reg::tensor2<double> out;

  out[0]          = mData[0];
  out[1] = out[2] = mData[1];
  out[3]          = mData[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2d<double> tensor2s<double>::cast<reg::tensor2d<double>>() const
{
  reg::tensor2d<double> out;

  out[0] = mData[0];
  out[1] = mData[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2<double> tensor2d<double>::cast<reg::tensor2<double>>() const
{
  reg::tensor2<double> out;

  out[0] = mData[0];
  out[3] = mData[1];

  out[1] = out[2] = 0.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2s<double> tensor2d<double>::cast<reg::tensor2s<double>>() const
{
  reg::tensor2s<double> out;

  out[0] = mData[0];
  out[2] = mData[1];

  out[1] = 0.0;

  return out;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X> inline size_t tensor4 <X>::size() const { return mSize; }
template<class X> inline size_t tensor2 <X>::size() const { return mSize; }
template<class X> inline size_t tensor2s<X>::size() const { return mSize; }
template<class X> inline size_t tensor2d<X>::size() const { return mSize; }
template<class X> inline size_t vector  <X>::size() const { return mSize; }
template<class X> inline size_t tensor4 <X>::ndim() const { return mNd;   }
template<class X> inline size_t tensor2 <X>::ndim() const { return mNd;   }
template<class X> inline size_t tensor2s<X>::ndim() const { return mNd;   }
template<class X> inline size_t tensor2d<X>::ndim() const { return mNd;   }
template<class X> inline size_t vector  <X>::ndim() const { return mNd;   }

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::shape() const
{
  std::vector<size_t> out(4);

  for ( size_t i = 0 ; i < 4 ; ++i )
    out[i] = mNd;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::shape() const
{
  std::vector<size_t> out(2);

  for ( size_t i = 0 ; i < 2 ; ++i )
    out[i] = mNd;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> out(1);

  for ( size_t i = 0 ; i < 1 ; ++i )
    out[i] = mNd;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { 8, 4, 2, 1 };

  if ( bytes )
  {
    out[0] *= sizeof(X);
    out[1] *= sizeof(X);
    out[2] *= sizeof(X);
    out[3] *= sizeof(X);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { 2, 1 };

  if ( bytes )
  {
    out[0] *= sizeof(X);
    out[1] *= sizeof(X);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { 1 };

  if ( bytes )
    out[0] *= sizeof(X);

  return out;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X>
inline const X& tensor4<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2d<X>::operator[](size_t i) const
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
inline const X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return mData[i*8+j*4+k*2+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  return mData[i*2+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator()(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  if ( i == 0 ) {
    if ( j == 0 ) return mData[0];
    else          return mData[1];
  }
  else {
    if ( j == 0 ) return mData[1];
    else          return mData[2];
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2d<X>::operator()(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  if (i == j) return mData[i];
  else        return mZero[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator()(size_t i) const
{
  assert( i < mNd );

  return mData[i];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline size_t tensor4<X>::compress(size_t i, size_t j, size_t k, size_t l) const
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return i*mNd*mNd*mNd+j*mNd*mNd+k*mNd+l;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2<X>::compress(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  return i*mNd+j;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2s<X>::compress(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  if (i <= j) return i*mNd - (i-1)*i/2 + j - i;
  else        return j*mNd - (j-1)*j/2 + i - j;
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X>
inline std::vector<size_t> tensor4<X>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(4);

  // reconstruct
  idx[0] = (i - i%(mNd*mNd*mNd)) / (mNd*mNd*mNd);  i -= idx[0] * (mNd*mNd*mNd);
  idx[1] = (i - i%(mNd*mNd)    ) / (mNd*mNd);      i -= idx[1] * (mNd*mNd);
  idx[2] = (i - i% mNd         ) /  mNd;           i -= idx[2] *  mNd;
  idx[3] =  i;

  return idx;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(2);

  // reconstruct
  idx[1] = i % mNd;
  idx[0] = ( i - idx[1] ) / mNd;

  return idx;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2s<X>::decompress(size_t i) const
{
  // check input
  assert( i < mSize );

  // allocate array-index
  std::vector<size_t> idx(2);

  // reconstruct
  idx[0] = 0;
  size_t keyafter;
  do
  {
    idx[0]++;
    keyafter = idx[0] * mNd - (idx[0] - 1) * idx[0] / 2;
  } while ( i >= keyafter );
  idx[0]--;
  idx[1] = mNd - keyafter + i;

  return idx;
}

// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X> inline const X* tensor4 <X>::data()  const { return std::begin(mData);         }
template<class X> inline auto     tensor4 <X>::begin() const { return std::begin(mData);         }
template<class X> inline auto     tensor4 <X>::end()   const { return std::begin(mData) + mSize; }
template<class X> inline const X* tensor2 <X>::data()  const { return std::begin(mData);         }
template<class X> inline auto     tensor2 <X>::begin() const { return std::begin(mData);         }
template<class X> inline auto     tensor2 <X>::end()   const { return std::begin(mData) + mSize; }
template<class X> inline const X* tensor2s<X>::data()  const { return std::begin(mData);         }
template<class X> inline auto     tensor2s<X>::begin() const { return std::begin(mData);         }
template<class X> inline auto     tensor2s<X>::end()   const { return std::begin(mData) + mSize; }
template<class X> inline const X* tensor2d<X>::data()  const { return std::begin(mData);         }
template<class X> inline auto     tensor2d<X>::begin() const { return std::begin(mData);         }
template<class X> inline auto     tensor2d<X>::end()   const { return std::begin(mData) + mSize; }
template<class X> inline const X* vector  <X>::data()  const { return std::begin(mData);         }
template<class X> inline auto     vector  <X>::begin() const { return std::begin(mData);         }
template<class X> inline auto     vector  <X>::end()   const { return std::begin(mData) + mSize; }

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto tensor4<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2s<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2d<X>::index(size_t i) const
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
inline auto tensor4<X>::item(size_t i, size_t j, size_t k, size_t l) const
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return begin() + i*mNd*mNd*mNd+j*mNd*mNd+k*mNd+l;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2<X>::item(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  return begin() + i*mNd+j;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2s<X>::item(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  if (i <= j) return begin() + i*mNd - (i-1)*i/2 + j - i;
  else        return begin() + j*mNd - (j-1)*j/2 + i - j;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::item(size_t i) const
{
  assert( i < mNd );

  return begin() + i;
}

// arithmetic operators
// =================================================================================================

template<class X>
inline reg::tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
    C[i+3] = A[i+3] * B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
    C[i+3] = A[i+3] / B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
    C[i+3] = A[i+3] + B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
    C[i+3] = A[i+3] - B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
    C[i+3] = A[i+3] * B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
    C[i+3] = A[i+3] / B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
    C[i+3] = A[i+3] + B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
    C[i+3] = A[i+3] - B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
    C[i+3] = A * B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
    C[i+3] = A / B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
    C[i+3] = A + B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
    C[i+3] = A - B[i+3];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[2] * B[1];
  C[3] = A[3] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[2] / B[1];
  C[3] = A[3] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[2] + B[1];
  C[3] = A[3] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  reg::tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[2] - B[1];
  C[3] = A[3] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  reg::tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[1] * B[2];
  C[3] = A[2] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[1] / B[2];
  C[3] = A[2] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[1] + B[2];
  C[3] = A[2] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  reg::tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[1] - B[2];
  C[3] = A[2] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  reg::tensor2<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] = A[1] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] = A[1] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] = A[1] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  reg::tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[2] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  reg::tensor2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator* (const vector<X> &A, const X &B)
{
  reg::vector<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator/ (const vector<X> &A, const X &B)
{
  reg::vector<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator+ (const vector<X> &A, const X &B)
{
  reg::vector<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator- (const vector<X> &A, const X &B)
{
  reg::vector<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator* (const X &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator/ (const X &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator+ (const X &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> operator- (const X &A, const vector<X> &B)
{
  reg::vector<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];

  return C;
}

// =================================================================================================
// tensor products: ddot
// =================================================================================================

template<class X>
inline reg::tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          for ( size_t m = 0 ; m < mNd ; ++m )
            for ( size_t n = 0 ; n < mNd ; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l) * B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  return mData[0]*B[0] + mData[3]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  X C;

  C  = mData[0] * B[0];
  C += mData[1] * B[1] * static_cast<X>(2);
  C += mData[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return mData[0]*B[0] + mData[2]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      for ( size_t l = 0 ; l < mNd ; ++l )
        C(k,l) += mData[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2<X> &B) const
{
  return mData[0]*B[0] + mData[1]*B[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  return mData[0]*B[0] + mData[1]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  return mData[0]*B[0] + mData[1]*B[1];
}


// =================================================================================================
// tensor products: dot
// =================================================================================================

template<class X>
inline reg::tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2d<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  reg::vector<X> C = reg::vector<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i) += (*this)(i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::dot(const vector<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// =================================================================================================
// tensor products: dyadic
// =================================================================================================

template<class X>
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      for ( size_t l = 0 ; l < mNd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      for ( size_t l = 0 ; l < mNd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C = reg::tensor4<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  reg::tensor2<X> C = reg::tensor2<X>::Zero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// =================================================================================================
// tensor products: ddot (out-of-class)
// =================================================================================================

template<class X>
inline reg::tensor4<X> ddot(const tensor4<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor4<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor4<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor4<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor2<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor2s<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> ddot(const tensor2d<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X ddot(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// =================================================================================================
// tensor products: dot (out-of-class)
// =================================================================================================

template<class X>
inline reg::tensor2<X> dot(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dot(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> dot(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const tensor2<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const tensor2s<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const tensor2d<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const vector<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const vector<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> dot(const vector<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X dot(const vector<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// =================================================================================================
// tensor products: dyadic (out-of-class)
// =================================================================================================

template<class X>
inline reg::tensor4<X> dyadic(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> dyadic(const vector<X> &A, const vector<X> &B)
{
  return A.dyadic(B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline reg::tensor4<X> tensor4<X>::T() const
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor4<X>::RT() const
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor4<X>::LT() const
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::T() const
{
  reg::tensor2<X> C;

  C[0] = mData[0];
  C[2] = mData[1];
  C[1] = mData[2];
  C[3] = mData[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> tensor2s<X>::T() const
{
  reg::tensor2s<X> C;

  C[0] = mData[0];
  C[1] = mData[1];
  C[2] = mData[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::T() const
{
  reg::tensor2d<X> C;

  C[0] = mData[0];
  C[1] = mData[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> transpose(const tensor2<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> transpose(const tensor2s<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> transpose(const tensor2d<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> transpose(const tensor4<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> transposeR(const tensor4<X> &A)
{
  return A.RT();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> transposeL(const tensor4<X> &A)
{
  return A.LT();
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline X tensor2<X>::trace() const
{
  return mData[0] + mData[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::trace() const
{
  return mData[0] + mData[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::trace() const
{
  return mData[0] + mData[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::det() const
{
  return mData[0] * mData[3] - mData[1] * mData[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::det() const
{
  return mData[0] * mData[2] - mData[1] * mData[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::det() const
{
  return mData[0] * mData[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  reg::tensor2<X> C;

  C[0] =                      mData[3] / D;
  C[1] = static_cast<X>(-1) * mData[1] / D;
  C[2] = static_cast<X>(-1) * mData[2] / D;
  C[3] =                      mData[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  reg::tensor2s<X> C;

  C[0] =                      mData[2] / D;
  C[1] = static_cast<X>(-1) * mData[1] / D;
  C[2] =                      mData[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  reg::tensor2d<X> C;

  C[0] = static_cast<X>(1) / mData[0];
  C[1] = static_cast<X>(1) / mData[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> inv(const tensor2<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> inv(const tensor2s<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> inv(const tensor2d<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X det(const tensor2<X> &A)
{
  return A.det();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X det(const tensor2s<X> &A)
{
  return A.det();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X det(const tensor2d<X> &A)
{
  return A.det();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X trace(const tensor2<X> &A)
{
  return A.trace();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X trace(const tensor2s<X> &A)
{
  return A.trace();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X trace(const tensor2d<X> &A)
{
  return A.trace();
}

// =================================================================================================
// equality operators
// =================================================================================================

template<class X>
inline bool tensor4<X>::operator== (const tensor4<X> &B) const
{
  for ( size_t i = 0 ; i < 16 ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor4<X>::operator== (const reg::tensor4<X> &B) const
{
  for ( size_t i = 0 ; i < 16 ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[2] ) return false;
  if ( mData[3] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[1] ) return false;
  if ( mData[3] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[3] != B[1] ) return false;
  if ( mData[1]         ) return false;
  if ( mData[2]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[2] ) return false;
  if ( mData[3] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const reg::tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[1] ) return false;
  if ( mData[3] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const reg::tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[3] != B[1] ) return false;
  if ( mData[1]         ) return false;
  if ( mData[2]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[1] != B[2] ) return false;
  if ( mData[2] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[2] != B[1] ) return false;
  if ( mData[1]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const reg::tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[2] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;
  if ( mData[1] != B[2] ) return false;
  if ( mData[2] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const reg::tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[2] != B[1] ) return false;
  if ( mData[1]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[3] ) return false;
  if (             B[1] ) return false;
  if (             B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[2] ) return false;
  if (             B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const reg::tensor2d<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[3] ) return false;
  if (             B[1] ) return false;
  if (             B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const reg::tensor2s<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[2] ) return false;
  if (             B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const vector<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const reg::vector<X> &B) const
{
  if ( mData[0] != B[0] ) return false;
  if ( mData[1] != B[1] ) return false;

  return true;
}

// =================================================================================================
// structure check
// =================================================================================================

template<class X>
inline bool tensor2<X>::issymmetric() const
{
  if ( mData[1] != mData[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::isdiagonal() const
{
  if ( mData[1] != static_cast<X>(0) ) return false;
  if ( mData[2] != static_cast<X>(0) ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::isdiagonal() const
{
  if ( mData[1] != static_cast<X>(0) ) return false;

  return true;
}

// =================================================================================================
// basic algebra: location of the minimum
// =================================================================================================

template<class X>
inline std::vector<size_t> tensor4<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2s<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2d<X>::argmin() const
{
  size_t i = std::min_element(begin(),end()) - begin();

  std::vector<size_t> out(mNd);

  std::fill(out.begin(), out.end(), i);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::argmin() const
{
  return decompress( std::min_element(begin(),end()) - begin() );
}

// =================================================================================================
// basic algebra: location of the minimum
// =================================================================================================

template<class X>
inline size_t tensor4<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2s<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2d<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::argminIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra: location of the maximum
// =================================================================================================

template<class X>
inline std::vector<size_t> tensor4<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2s<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2d<X>::argmax() const
{
  size_t i = std::max_element(begin(),end()) - begin();

  std::vector<size_t> out(mNd);

  std::fill(out.begin(), out.end(), i);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::argmax() const
{
  return decompress( std::max_element(begin(),end()) - begin() );
}

// =================================================================================================
// basic algebra: location of the maximum
// =================================================================================================

template<class X>
inline size_t tensor4<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2s<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2d<X>::argmaxIndex() const
{
  return std::max_element(begin(),end()) - begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::argmaxIndex() const
{
  return std::min_element(begin(),end()) - begin();
}

// =================================================================================================
// basic algebra: minimum
// =================================================================================================

template<class X>
inline X tensor4<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::minCoeff() const
{
  return *std::min_element(begin(),end());
}

// =================================================================================================
// basic algebra: maximum
// =================================================================================================

template<class X>
inline X tensor4<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::maxCoeff() const
{
  return *std::max_element(begin(),end());
}

// =================================================================================================
// basic algebra: sum
// =================================================================================================

template<class X>
inline X tensor4<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i];

  return out;
}

// =================================================================================================
// basic algebra: mean
// =================================================================================================

template<class X>
inline double tensor4<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2s<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2d<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double vector<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(mSize);
}

// =================================================================================================
// basic algebra: weighted average
// =================================================================================================

template<class X>
inline double tensor4<X>::average(const tensor4<X> &weights, bool norm) const
{
  assert( mNd = weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2<X>::average(const tensor2<X> &weights, bool norm) const
{
  assert( mNd = weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2s<X>::average(const tensor2s<X> &weights, bool norm) const
{
  assert( mNd = weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double tensor2d<X>::average(const tensor2d<X> &weights, bool norm) const
{
  assert( mNd = weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double vector<X>::average(const vector<X> &weights, bool norm) const
{
  assert( mNd = weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    out += mData[i] * weights[i];

  if ( norm ) return static_cast<double>(out)/static_cast<double>(weights.sum());
  else        return static_cast<double>(out);
}

// =================================================================================================
// basic algebra: norm
// =================================================================================================

template<class X>
inline X tensor4<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C += std::abs(mData[i  ]);
    C += std::abs(mData[i+1]);
    C += std::abs(mData[i+2]);
    C += std::abs(mData[i+3]);
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::norm() const
{
  X C;

  C  = std::abs(mData[0]);
  C += std::abs(mData[1]);
  C += std::abs(mData[2]);
  C += std::abs(mData[3]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::norm() const
{
  X C;

  C  = std::abs(mData[0]);
  C += std::abs(mData[1]);
  C += std::abs(mData[2]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::norm() const
{
  X C;

  C  = std::abs(mData[0]);
  C += std::abs(mData[1]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::norm() const
{
  X C;

  C  = std::abs(mData[0]);
  C += std::abs(mData[1]);

  return C;
}

// =================================================================================================
// basic algebra: length
// =================================================================================================

template<class X>
inline X vector<X>::length() const
{
  X C;

  C  = std::pow(mData[0],2.);
  C += std::pow(mData[1],2.);

  return std::sqrt(C);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void tensor4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(mNd).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, (*this)(i,j,k,l) );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(0,0), (*this)(0,1));
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(1,0), (*this)(1,1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(0,0), (*this)(0,1));
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(1,0), (*this)(1,1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(0,0), (*this)(0,1));
  std::printf((fmt + "," + fmt + ";\n").c_str(), (*this)(1,0), (*this)(1,1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + "\n").c_str(), mData[0], mData[1]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor4<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < 2 ; ++i ) {
    for ( size_t j = 0 ; j < 2 ; ++j ) {
      for ( size_t k = 0 ; k < 2 ; ++k ) {
        for ( size_t l = 0 ; l < 2 ; ++l ) {
          out << "(" << i << "," << j << "," << k << "," << l << ") = ";
          out << std::setw(w) << std::setprecision(p) << src(i,j,k,l);
          if ( !(i==1 and j==1 and k==1 and l==1) ) out << std::endl;
        }
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor2<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor2s<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor2d<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1);

  return out;
}

// =================================================================================================

}}} // namespace ...

#endif

