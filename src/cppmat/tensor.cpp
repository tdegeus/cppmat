/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR_CPP
#define CPPMAT_TENSOR_CPP

// -------------------------------------------------------------------------------------------------

#include "tensor.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline tensor4<X>::tensor4(size_t nd)
{
  // change size of internal storage
  resize(nd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Arange(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Zero(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Ones(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Constant(size_t nd, X D)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::I(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Irt(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setIrt();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Is(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setIs();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Id(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setId();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::II(size_t nd)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setII();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor4<X> tensor4<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor4<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>::tensor2(size_t nd)
{
  // change size of internal storage
  resize(nd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Arange(size_t nd)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Zero(size_t nd)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Ones(size_t nd)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Constant(size_t nd, X D)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::I(size_t nd)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2<X> tensor2<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor2<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>::tensor2s(size_t nd)
{
  // change size of internal storage
  resize(nd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Arange(size_t nd)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Zero(size_t nd)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Ones(size_t nd)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Constant(size_t nd, X D)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::I(size_t nd)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2s<X> tensor2s<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2s<X> tensor2s<X>::CopyDense(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor2s<X> out(nd);

  // initialize
  out.setCopyDense(first,last);

  return out;
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
inline tensor2d<X>::tensor2d(size_t nd)
{
  // change size of internal storage
  resize(nd);

  // set dummy parameter
  mZero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Arange(size_t nd)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Zero(size_t nd)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Ones(size_t nd)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Constant(size_t nd, X D)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::I(size_t nd)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2d<X> tensor2d<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2d<X> tensor2d<X>::CopyDense(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  tensor2d<X> out(nd);

  // initialize
  out.setCopyDense(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector(size_t nd)
{
  // change size of internal storage
  resize(nd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Arange(size_t nd)
{
  // call basic constructor
  vector<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Zero(size_t nd)
{
  // call basic constructor
  vector<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Ones(size_t nd)
{
  // call basic constructor
  vector<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Constant(size_t nd, X D)
{
  // call basic constructor
  vector<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline vector<X> vector<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  // call basic constructor
  vector<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// cast to different class
// =================================================================================================

template<>
template<>
inline tensor2s<double> tensor2<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i ; j < mNd ; ++j )
      out[i*mNd-(i-1)*i/2+j-i] = ( mData[i*mNd+j] + mData[j*mNd+i] ) / 2.;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i] = mData[i*mNd+i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2s<double>::cast<tensor2<double>>() const
{
  tensor2<double> out(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i ; j < mNd ; ++j )
      out[i*mNd+j] = out[j*mNd+i] = mData[i*mNd-(i-1)*i/2+j-i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2s<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i] = mData[i*mNd-(i-1)*i/2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2d<double>::cast<tensor2<double>>() const
{
  tensor2<double> out = tensor2<double>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i*mNd+i] = mData[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2s<double> tensor2d<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out = tensor2s<double>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i*mNd-(i-1)*i/2] = mData[i];

  return out;
}

// =================================================================================================
// automatic cast to different class
// =================================================================================================

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2s<X>::operator tensor2<X> () const
{
  tensor2<X> out(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i ; j < mNd ; ++j )
      out[i*mNd+j] = out[j*mNd+i] = mData[i*mNd-(i-1)*i/2+j-i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2<X> () const
{
  tensor2<double> out = tensor2<double>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i*mNd+i] = mData[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2s<X> () const
{
  tensor2s<X> out(mNd,static_cast<X>(0));

  for ( size_t i = 0 ; i < mNd ; ++i )
    out[i*mNd-(i-1)*i/2] = mData[i];

  return out;
}
#endif

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void tensor4<X>::resize(size_t nd)
{
  // store size
  mNd    = nd;
  mSize = mNd*mNd*mNd*mNd;

  // resize storage
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::resize(size_t nd)
{
  // store size
  mNd    = nd;
  mSize = mNd*mNd;

  // resize storage
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::resize(size_t nd)
{
  // store size
  mNd    = nd;
  mSize = (mNd+1)*mNd/2;

  // resize storage
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::resize(size_t nd)
{
  // store size
  mNd    = nd;
  mSize = mNd;

  // resize storage
  mData.resize(mSize);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::resize(size_t nd)
{
  // store size
  mNd    = nd;
  mSize = mNd;

  // resize storage
  mData.resize(mSize);
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

  std::fill(out.begin(), out.end(), mNd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::shape() const
{
  std::vector<size_t> out(2);

  std::fill(out.begin(), out.end(), mNd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> out(1);

  std::fill(out.begin(), out.end(), mNd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { mNd*mNd*mNd, mNd*mNd, mNd, 1 };

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
  std::vector<size_t> out = { mNd, 1 };

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
inline X& tensor4<X>::operator[](size_t i)
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor4<X>::operator[](size_t i) const
{
  assert( i < mSize );

  return mData[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2<X>::operator[](size_t i)
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
inline X& tensor2s<X>::operator[](size_t i)
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
inline X& tensor2d<X>::operator[](size_t i)
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
inline X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l)
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return mData[i*mNd*mNd*mNd+j*mNd*mNd+k*mNd+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return mData[i*mNd*mNd*mNd+j*mNd*mNd+k*mNd+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2<X>::operator()(size_t i, size_t j)
{
  assert( i < mNd );
  assert( j < mNd );

  return mData[i*mNd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  return mData[i*mNd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2s<X>::operator()(size_t i, size_t j)
{
  assert( i < mNd );
  assert( j < mNd );

  if (i <= j) return mData[ i*mNd - (i-1)*i/2 + j - i ];
  else        return mData[ j*mNd - (j-1)*j/2 + i - j ];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator()(size_t i, size_t j) const
{
  assert( i < mNd );
  assert( j < mNd );

  if (i <= j) return mData[ i*mNd - (i-1)*i/2 + j - i ];
  else        return mData[ j*mNd - (j-1)*j/2 + i - j ];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2d<X>::operator()(size_t i, size_t j)
{
  assert( i < mNd );
  assert( j < mNd );

  if (i == j) return mData[i];
  else        return mZero[0];
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
inline X& vector<X>::operator()(size_t i)
{
  assert( i < mNd );

  return mData[i];
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

template<class X> inline X*       tensor4 <X>::data()        { return mData.data();  }
template<class X> inline const X* tensor4 <X>::data()  const { return mData.data();  }
template<class X> inline auto     tensor4 <X>::begin()       { return mData.begin(); }
template<class X> inline auto     tensor4 <X>::begin() const { return mData.begin(); }
template<class X> inline auto     tensor4 <X>::end()         { return mData.end();   }
template<class X> inline auto     tensor4 <X>::end()   const { return mData.end();   }
template<class X> inline X*       tensor2 <X>::data()        { return mData.data();  }
template<class X> inline const X* tensor2 <X>::data()  const { return mData.data();  }
template<class X> inline auto     tensor2 <X>::begin()       { return mData.begin(); }
template<class X> inline auto     tensor2 <X>::begin() const { return mData.begin(); }
template<class X> inline auto     tensor2 <X>::end()         { return mData.end();   }
template<class X> inline auto     tensor2 <X>::end()   const { return mData.end();   }
template<class X> inline X*       tensor2s<X>::data()        { return mData.data();  }
template<class X> inline const X* tensor2s<X>::data()  const { return mData.data();  }
template<class X> inline auto     tensor2s<X>::begin()       { return mData.begin(); }
template<class X> inline auto     tensor2s<X>::begin() const { return mData.begin(); }
template<class X> inline auto     tensor2s<X>::end()         { return mData.end();   }
template<class X> inline auto     tensor2s<X>::end()   const { return mData.end();   }
template<class X> inline X*       tensor2d<X>::data()        { return mData.data();  }
template<class X> inline const X* tensor2d<X>::data()  const { return mData.data();  }
template<class X> inline auto     tensor2d<X>::begin()       { return mData.begin(); }
template<class X> inline auto     tensor2d<X>::begin() const { return mData.begin(); }
template<class X> inline auto     tensor2d<X>::end()         { return mData.end();   }
template<class X> inline auto     tensor2d<X>::end()   const { return mData.end();   }
template<class X> inline X*       vector  <X>::data()        { return mData.data();  }
template<class X> inline const X* vector  <X>::data()  const { return mData.data();  }
template<class X> inline auto     vector  <X>::begin()       { return mData.begin(); }
template<class X> inline auto     vector  <X>::begin() const { return mData.begin(); }
template<class X> inline auto     vector  <X>::end()         { return mData.end();   }
template<class X> inline auto     vector  <X>::end()   const { return mData.end();   }

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto tensor4<X>::index(size_t i)
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor4<X>::index(size_t i) const
{
  assert( i < mSize );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto tensor2<X>::index(size_t i)
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
inline auto tensor2s<X>::index(size_t i)
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
inline auto tensor2d<X>::index(size_t i)
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
inline auto tensor4<X>::item(size_t i, size_t j, size_t k, size_t l)
{
  assert( i < mNd );
  assert( j < mNd );
  assert( k < mNd );
  assert( l < mNd );

  return begin() + i*mNd*mNd*mNd+j*mNd*mNd+k*mNd+l;
}

// -------------------------------------------------------------------------------------------------

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
inline auto tensor2<X>::item(size_t i, size_t j)
{
  assert( i < mNd );
  assert( j < mNd );

  return begin() + i*mNd+j;
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
inline auto tensor2s<X>::item(size_t i, size_t j)
{
  assert( i < mNd );
  assert( j < mNd );

  if (i <= j) return begin() + i*mNd - (i-1)*i/2 + j - i;
  else        return begin() + j*mNd - (j-1)*j/2 + i - j;
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
inline auto vector<X>::item(size_t i)
{
  assert( i < mNd );

  return begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto vector<X>::item(size_t i) const
{
  assert( i < mNd );

  return begin() + i;
}

// =================================================================================================
// basic initialization - tensor4
// =================================================================================================

template<class X>
inline void tensor4<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor4<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setI()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          if ( i == l and j == k )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIrt()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          if ( i == k and j == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIs()
{
  return ( tensor4<X>::I(mNd) + tensor4<X>::Irt(mNd) ) / static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setId()
{
  return tensor4<X>::Is(mNd) - tensor4<X>::II(mNd)/static_cast<X>(mNd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setII()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          if ( i == j and k == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// =================================================================================================
// basic initialization - tensor2
// =================================================================================================

template<class X>
inline void tensor2<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setI()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// basic initialization - tensor2s
// =================================================================================================

template<class X>
inline void tensor2s<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( mNd*mNd == last-first );

  // check for symmetry
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i+1 ; j < mNd ; ++j )
      assert( first[i*mNd+j] == first[j*mNd+i] );
  #endif

  // copy from input (ignores lower diagonal terms)
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i ; j < mNd ; ++j )
      mData[i*mNd-(i-1)*i/2+j-i] = first[i*mNd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setI()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// basic initialization - tensor2d
// =================================================================================================

template<class X>
inline void tensor2d<X>::setArange()
{
  for ( size_t i = 0 ; i < mSize ; ++i ) mData[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setZero()
{
  std::fill(begin(), end(), static_cast<X>(0));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setOnes()
{
  std::fill(begin(), end(), static_cast<X>(1));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setConstant(X D)
{
  std::fill(begin(), end(), D);
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopy(Iterator first, Iterator last)
{
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( mNd*mNd == last-first );

  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      if ( i !=j )
        assert( !first[i*mNd+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] = first[i*mNd+i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setI()
{
  setZero();

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] = static_cast<X>(1);
}

// =================================================================================================
// basic initialization - vector
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
  assert( mSize == last-first );

  std::copy(first, last, begin());
}

// =================================================================================================
// arithmetic operators - tensor4
// =================================================================================================

template<class X>
inline tensor4<X>& tensor4<X>::operator*= (const tensor4<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const tensor4<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const tensor4<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const tensor4<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2
// =================================================================================================

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      // - extract value
      X b = B[i*mNd-(i-1)*i/2+j-i];
      // - store symmetrically
                    mData[i*mNd+j] *= b;
      if ( i != j ) mData[j*mNd+i] *= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      // - extract value
      X b = B[i*mNd-(i-1)*i/2+j-i];
      // - store symmetrically
                    mData[i*mNd+j] /= b;
      if ( i != j ) mData[j*mNd+i] /= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      // - extract value
      X b = B[i*mNd-(i-1)*i/2+j-i];
      // - store symmetrically
                    mData[i*mNd+j] += b;
      if ( i != j ) mData[j*mNd+i] += b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      // - extract value
      X b = B[i*mNd-(i-1)*i/2+j-i];
      // - store symmetrically
                    mData[i*mNd+j] -= b;
      if ( i != j ) mData[j*mNd+i] -= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = 0 ; j < mNd ; ++j ) {
      if ( i == j ) mData[i*mNd+i] *= B[i];
      else          mData[i*mNd+j]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i*mNd+i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i*mNd+i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2s
// =================================================================================================

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j=i; j<mNd; ++j ) {
      if ( i == j ) mData[i*mNd-(i-1)*i/2    ] *= B[i];
      else          mData[i*mNd-(i-1)*i/2+j-i]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i*mNd-(i-1)*i/2] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i*mNd-(i-1)*i/2] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2d
// =================================================================================================

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const tensor2d<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd; ++i )
    mData[i] *= B[i*mNd+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd; ++i )
    mData[i] /= B[i*mNd+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] *= B[i*mNd-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2s<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] /= B[i*mNd-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    mData[i] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - vector
// =================================================================================================

template<class X>
inline vector<X>& vector<X>::operator*= (const vector<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const vector<X> &B)
{
  assert( mNd == B.ndim() );

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

// =================================================================================================
// arithmetic operators - mixed class
// =================================================================================================

template<class X>
inline tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  tensor4<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  tensor4<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  assert( A.size() == B.size() );
  assert( B.ndim() == B.ndim() );

  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] * b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] * b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] / b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] / b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] + b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] + b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X b = B[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = A[ i*nd + j ] - b;
      if ( i != j ) C[ j*nd + i ] = A[ j*nd + i ] - b;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i*nd + j ] + B[ i ];
      else          C[ i*nd + j ] = A[ i*nd + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i*nd + j ] - B[ i ];
      else          C[ i*nd + j ] = A[ i*nd + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  tensor2<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a * B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a * B[ j*nd + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a / B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a / B[ j*nd + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a + B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a + B[ j*nd + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      // - extract value
      X a = A[ i*nd - (i-1)*i/2 + j - i ];
      // - perform multiplication
                    C[ i*nd + j ] = a - B[ i*nd + j ];
      if ( i != j ) C[ j*nd + i ] = a - B[ j*nd + i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i ] + B[ i*nd + j ];
      else          C[ i*nd + j ] =          B[ i*nd + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd + j ] = A[ i ] - B[ i*nd + j ];
      else          C[ i*nd + j ] =        - B[ i*nd + j ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  tensor2<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i*nd - (i-1)*i/2         ] + B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i*nd - (i-1)*i/2         ] - B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] + B;
      else          C[ i*nd - (i-1)*i/2 + j - i ] =          B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] - B;
      else          C[ i*nd - (i-1)*i/2 + j - i ] =        - B;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] + B[ i*nd - (i-1)*i/2         ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] =          B[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A[ i ] - B[ i*nd - (i-1)*i/2         ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] =        - B[ i*nd - (i-1)*i/2 + j - i ];
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  size_t      nd = B.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A + B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  size_t      nd = B.ndim();
  tensor2s<X> C(nd);

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = i ; j < nd ; ++j ) {
      if ( i == j ) C[ i*nd - (i-1)*i/2         ] = A - B[ i ];
      else          C[ i*nd - (i-1)*i/2 + j - i ] = A;
    }
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[i] * B[ i*nd + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[i] / B[ i*nd + i ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[i] * B[ i*nd - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[i] / B[ i*nd - (i-1)*i/2 ];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[ i*nd + i ] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  assert( A.ndim() == B.ndim() );

  size_t nd = A.ndim();
  tensor2d<X> C(A.ndim());

  for ( size_t i = 0 ; i < nd ; ++i )
    C[i] = A[ i*nd - (i-1)*i/2 ] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  tensor2d<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.ndim() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C(A.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C(B.ndim());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

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
inline tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

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
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

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
  assert( ndim() == B.ndim() );

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
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

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
  assert( ndim() == B.ndim() );

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
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += ( mData[i*mNd-(i-1)*i/2] * B[i*mNd-(i-1)*i/2] );

  for ( size_t i = 0 ; i<mNd ; ++i )
    for ( size_t j = i+1 ; j<mNd ; ++j )
      C += ( static_cast<X>(2) * mData[i*mNd-(i-1)*i/2+j-i] * B[i*mNd-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += mData[i*mNd-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

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
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += mData[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += mData[i]*B[i*mNd-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += mData[i]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2d<X> C = tensor2d<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C(i) += (*this)(i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::dot(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      for ( size_t l = 0 ; l < mNd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      for ( size_t l = 0 ; l < mNd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t k = 0 ; k < mNd ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::cross(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  if ( mNd != 3 )
    throw std::runtime_error("'cross' only implemented in 3D");

  vector<X> C(3);

  C[0] =                     mData[1]*B[2]-B[1]*mData[2] ;
  C[1] = static_cast<X>(-1)*(mData[0]*B[2]-B[0]*mData[2]);
  C[2] =                     mData[0]*B[1]-B[0]*mData[1] ;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> ddot(const tensor4<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor4<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor4<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor4<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor2<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor2s<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> ddot(const tensor2d<X> &A, const tensor4<X> &B)
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

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dot(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> dot(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const tensor2<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const tensor2s<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const tensor2d<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const vector<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const vector<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> dot(const vector<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X dot(const vector<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> dyadic(const vector<X> &A, const vector<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> cross (const vector<X> &A, const vector<X> &B)
{
  return A.cross (B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::T() const
{
  tensor4<X> C(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::RT() const
{
  tensor4<X> C(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::LT() const
{
  tensor4<X> C(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      for ( size_t k = 0 ; k < mNd ; ++k )
        for ( size_t l = 0 ; l < mNd ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::T() const
{
  tensor2<X> C(mNd);

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      C(j,i) = (*this)(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::T() const
{
  tensor2s<X> C(mNd);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C[i] = (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::T() const
{
  tensor2d<X> C(mNd);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C[i] = (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> transpose(const tensor2<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> transpose(const tensor2s<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> transpose(const tensor2d<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> transpose(const tensor4<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> transposeR(const tensor4<X> &A)
{
  return A.RT();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> transposeL(const tensor4<X> &A)
{
  return A.LT();
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline X tensor2<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::det() const
{
  if ( mNd==2 )
   return mData[0] * mData[3] - mData[1] * mData[2];

  if ( mNd==3 )
    return ( mData[0] * mData[4] * mData[8] +
             mData[1] * mData[5] * mData[6] +
             mData[2] * mData[3] * mData[7] ) -
           ( mData[2] * mData[4] * mData[6] +
             mData[1] * mData[3] * mData[8] +
             mData[0] * mData[5] * mData[7] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::det() const
{
  if ( mNd==2 )
   return mData[0] * mData[2] - mData[1] * mData[1];

  if ( mNd==3 )
    return (                     mData[0] * mData[3] * mData[5] +
             static_cast<X>(2) * mData[1] * mData[2] * mData[4] ) -
           (                     mData[4] * mData[4] * mData[0] +
                                 mData[2] * mData[2] * mData[3] +
                                 mData[1] * mData[1] * mData[5] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::det() const
{
  X C = static_cast<X>(1);

  for ( size_t i = 0 ; i < mNd ; ++i )
    C *= (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2<X> C(mNd);

  if ( mNd==2 )
  {
    C[0] =                      mData[3] / D;
    C[1] = static_cast<X>(-1) * mData[1] / D;
    C[2] = static_cast<X>(-1) * mData[2] / D;
    C[3] =                      mData[0] / D;

    return C;
  }

  if ( mNd==3 )
  {
    C[0] = (mData[4]*mData[8]-mData[5]*mData[7]) / D;
    C[1] = (mData[2]*mData[7]-mData[1]*mData[8]) / D;
    C[2] = (mData[1]*mData[5]-mData[2]*mData[4]) / D;
    C[3] = (mData[5]*mData[6]-mData[3]*mData[8]) / D;
    C[4] = (mData[0]*mData[8]-mData[2]*mData[6]) / D;
    C[5] = (mData[2]*mData[3]-mData[0]*mData[5]) / D;
    C[6] = (mData[3]*mData[7]-mData[4]*mData[6]) / D;
    C[7] = (mData[1]*mData[6]-mData[0]*mData[7]) / D;
    C[8] = (mData[0]*mData[4]-mData[1]*mData[3]) / D;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2s<X> C(mNd);

  if ( mNd==2 )
  {
    C[0] =                      mData[2] / D;
    C[1] = static_cast<X>(-1) * mData[1] / D;
    C[2] =                      mData[0] / D;

    return C;
  }

  if ( mNd==3 )
  {
    C[0] = (mData[3]*mData[5]-mData[4]*mData[4]) / D;
    C[1] = (mData[2]*mData[4]-mData[1]*mData[5]) / D;
    C[2] = (mData[1]*mData[4]-mData[2]*mData[3]) / D;
    C[3] = (mData[0]*mData[5]-mData[2]*mData[2]) / D;
    C[4] = (mData[2]*mData[1]-mData[0]*mData[4]) / D;
    C[5] = (mData[0]*mData[3]-mData[1]*mData[1]) / D;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C(mNd);

  for ( size_t i = 0; i < mNd ; ++i )
    C[i] = static_cast<X>(1) / mData[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> inv(const tensor2<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> inv(const tensor2s<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> inv(const tensor2d<X> &A)
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
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2s<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      if ( mData[i*mNd+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2d<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      if ( mData[i*mNd+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      if ( mData[i*mNd-(i-1)*i/2+j-i] != B(i,j) ) return false;
      if ( mData[i*mNd-(i-1)*i/2+j-i] != B(j,i) ) return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i ; j < mNd ; ++j )
      if ( mData[i*mNd-(i-1)*i/2+j-i] != B(i,j) ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2d<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = 0 ; j < mNd ; ++j ) {
      if ( i == j ) { if ( mData[i] != B(i,i) ) return false; }
      else          { if (              B(i,j) ) return false; }
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2s<X> &B) const
{
  assert( mNd == B.ndim() );

  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = i ; j < mNd ; ++j ) {
      if ( i == j ) { if ( mData[i] != B(i,i) ) return false; }
      else          { if (              B(i,j) ) return false; }
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const vector<X> &B) const
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    if ( mData[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================
// structure check
// =================================================================================================

template<class X>
inline bool tensor2<X>::issymmetric() const
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i+1 ; j < mNd ; ++j )
      if ( mData[ i*mNd + j ] != mData[ j*mNd + i ] )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = 0 ; j < mNd ; ++j )
      if ( i != j )
        if ( mData[ i*mNd + j ] )
          return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < mNd ; ++i )
    for ( size_t j = i+1 ; j < mNd ; ++j )
      if ( mData[i*mNd-(i-1)*i/2+j-i] )
        return false;

  return true;
}

// =================================================================================================
// basic algebra: absolute value
// =================================================================================================

template<class X>
inline void tensor4<X>::abs()
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = std::abs(mData[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::abs()
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = std::abs(mData[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::abs()
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = std::abs(mData[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::abs()
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = std::abs(mData[i]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::abs()
{
  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] = std::abs(mData[i]);
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

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::abs(mData[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::abs(mData[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::abs(mData[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::abs(mData[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::abs(mData[i]);

  return C;
}

// =================================================================================================
// basic algebra: length
// =================================================================================================

template<class X>
inline X vector<X>::length() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < mSize ; ++i )
    C += std::pow(mData[i],2.);

  return std::sqrt(C);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setUnitLength()
{
  X C = length();

  if ( C <= static_cast<X>(0) ) return;

  for ( size_t i = 0 ; i < mSize ; ++i )
    mData[i] /= C;
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
  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = 0 ; j < mNd ; ++j ) {
      if ( j != mNd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else              std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = 0 ; j < mNd ; ++j ) {
      if ( j != mNd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else              std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < mNd ; ++i ) {
    for ( size_t j = 0 ; j < mNd ; ++j ) {
      if ( j != mNd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else              std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < mNd ; ++j ) {
    if ( j != mNd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else              std::printf((fmt + ";\n").c_str(), (*this)(j));
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor4<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  size_t nd = src.ndim();

  for ( size_t i = 0 ; i < nd ; ++i ) {
    for ( size_t j = 0 ; j < nd ; ++j ) {
      for ( size_t k = 0 ; k < nd ; ++k ) {
        for ( size_t l = 0 ; l < nd ; ++l ) {
          out << "(" << i << "," << j << "," << k << "," << l << ") = ";
          out << std::setw(w) << std::setprecision(p) << src(i,j,k,l);
          if ( !(i==nd-1 and j==nd-1 and k==nd-1 and l==nd-1) ) out << std::endl;
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

  for ( size_t i = 0 ; i < src.ndim() ; ++i ) {
    for ( size_t j = 0 ; j < src.ndim() ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != src.ndim()-1 ) out << ", ";
      else if ( i != src.ndim()-1 ) out << ";" << std::endl;
      else                          out << ";";
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor2s<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < src.ndim() ; ++i ) {
    for ( size_t j = 0 ; j < src.ndim() ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != src.ndim()-1 ) out << ", ";
      else if ( i != src.ndim()-1 ) out << ";" << std::endl;
      else                          out << ";";
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const tensor2d<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < src.ndim() ; ++i ) {
    for ( size_t j = 0 ; j < src.ndim() ; ++j ) {
      out << std::setw(w) << std::setprecision(p) << src(i,j);
      if      ( j != src.ndim()-1 ) out << ", ";
      else if ( i != src.ndim()-1 ) out << ";" << std::endl;
      else                          out << ";";
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, const vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t j = 0 ; j < src.ndim() ; ++j ) {
    out << std::setw(w) << std::setprecision(p) << src(j);
    if ( j != src.ndim()-1 ) out << ", ";
  }

  return out;
}

// =================================================================================================

}} // namespace ...

#endif

