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
  // store shape, and other size parameters, allocate "m_data"
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
  m_ndim = shape.size();

  // initialize shape in all directions to "1"
  for ( size_t i = 0 ; i < MAX_DIM ; ++i )
  {
    m_shape  [i] = 1;
    m_strides[i] = 1;
  }

  // initialize size
  m_size = 1;

  // update shape/size
  for ( size_t i = 0 ; i < m_ndim ; ++i )
  {
    m_shape[i] = shape[i];
    m_size    *= shape[i];
  }

  // set storage strides
  for ( size_t i = 0 ; i < m_ndim ; ++i )
    for ( size_t j = i+1 ; j < m_ndim ; ++j )
      m_strides[i] *= m_shape[j];

  // allocate data
  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::chdim(size_t ndim)
{
  // check that all removed dimensions are of shape 1
  #ifndef NDEBUG
    for ( size_t i = ndim ; i < MAX_DIM ; ++i ) assert( m_shape[i] == 1 );
  #endif

  // update number of dimensions
  m_ndim = ndim;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::reshape(const std::vector<size_t> &shape)
{
  // check that the size is unchanged
  #ifndef NDEBUG
    size_t n = 1;

    for ( auto &i : shape ) n *= i;

    assert( n == m_size );
  #endif

  // process new shape
  resize(shape);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t array<X>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::shape(int i) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( i  <      static_cast<int>(m_ndim) );
  assert( i >= -1 * static_cast<int>(m_ndim) );

  // get number of dimensions as integer
  int n = static_cast<int>(m_ndim);

  // correct periodic index
  i = ( n + (i%n) ) % n;

  // return shape
  return m_shape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::shape(size_t i) const
{
  // check axis: (0,1,...,ndim-1)
  assert( i < m_ndim );

  // return shape
  return m_shape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> array<X>::shape() const
{
  std::vector<size_t> ret(m_ndim);

  for ( size_t i = 0 ; i < m_ndim ; ++i ) ret[i] = m_shape[i];

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> array<X>::strides(bool bytes) const
{
  std::vector<size_t> ret(m_ndim);

  for ( size_t i = 0 ; i < m_ndim ; ++i )
    ret[i] = m_strides[i];

  if ( bytes )
    for ( size_t i = 0 ; i < m_ndim ; ++i )
      ret[i] *= sizeof(X);

  return ret;
}

// =================================================================================================
// index operators : operator[...]
// =================================================================================================

template<class X>
inline X& array<X>::operator[](size_t i)
{
  assert( i < m_size );

  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator[](size_t i) const
{
  assert( i < m_size );

  return m_data[i];
}

// =================================================================================================
// index operators : operator(...)
// =================================================================================================

template<class X>
inline X& array<X>::operator()(size_t a)
{
  assert( m_ndim >= 1 );

  assert( a < m_shape[0] );

  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a) const
{
  assert( m_ndim >= 1 );

  assert( a < m_shape[0] );

  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b)
{
  assert( m_ndim >= 2 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );

  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b) const
{
  assert( m_ndim >= 2 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );

  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c)
{
  assert( m_ndim >= 3 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c) const
{
  assert( m_ndim >= 3 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d)
{
  assert( m_ndim >= 4 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d) const
{
  assert( m_ndim >= 4 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( m_ndim >= 5 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( m_ndim >= 5 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& array<X>::operator()(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( m_ndim >= 6 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );
  assert( f < m_shape[5] );

  return \
  m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& array<X>::operator()(
  size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( m_ndim >= 6 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );
  assert( f < m_shape[5] );

  return \
  m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// =================================================================================================
// index operators : at(...)
// =================================================================================================

template<class X>
template<class Iterator>
inline X& array<X>::at(Iterator first, Iterator last)
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0      );
  assert( static_cast<size_t>(last-first) <= m_ndim );

  // iterator to shape and stride
  size_t *shape  = &m_shape  [0];
  size_t *stride = &m_strides[0];

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

  return m_data[idx];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline const X& array<X>::at(Iterator first, Iterator last) const
{
  // check input
  assert( static_cast<size_t>(last-first)  > 0      );
  assert( static_cast<size_t>(last-first) <= m_ndim );

  // iterator to shape and stride
  size_t *shape  = &m_shape  [0];
  size_t *stride = &m_strides[0];

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

  return m_data[idx];
}

// =================================================================================================
// index operators : compress(...)
// =================================================================================================

template<class X>
inline size_t array<X>::compress(size_t a) const
{
  assert( m_ndim >= 1 );

  assert( a < m_shape[0] );

  return a*m_strides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b) const
{
  assert( m_ndim >= 2 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );

  return a*m_strides[0]+b*m_strides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c) const
{
  assert( m_ndim >= 3 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );

  return a*m_strides[0]+b*m_strides[1]+c*m_strides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d) const
{
  assert( m_ndim >= 4 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );

  return a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( m_ndim >= 5 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );

  return a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t array<X>::compress(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( m_ndim >= 6 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );
  assert( f < m_shape[5] );

  return a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5];
}

// =================================================================================================
// index operators : decompress(...)
// =================================================================================================

template<class X>
inline std::vector<size_t> array<X>::decompress(size_t i) const
{
  // check input
  assert( i < m_size );

  // allocate array-index
  std::vector<size_t> idx(m_ndim);

  // reconstruct
  for ( size_t j = 0 ; j < m_ndim ; ++j ) {
    idx[j] = (i - i%m_strides[j]) / m_strides[j];
    i -= idx[j] * m_strides[j];
  }

  return idx;
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline X* array<X>::data()
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* array<X>::data() const
{
  return m_data.data();
}

// =================================================================================================
// iterators : begin() and end()
// =================================================================================================

template<class X>
inline auto array<X>::begin()
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::end()
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::end() const
{
  return m_data.end();
}

// =================================================================================================
// iterators : index()
// =================================================================================================

template<class X>
inline auto array<X>::index(size_t i)
{
  assert( i < m_size );

  return m_data.begin() + i;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::index(size_t i) const
{
  assert( i < m_size );

  return m_data.begin() + i;
}

// =================================================================================================
// iterators : item()
// =================================================================================================

template<class X>
inline auto array<X>::item(size_t a)
{
  assert( m_ndim >= 1 );

  assert( a < m_shape[0] );

  return m_data.begin() + a*m_strides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a) const
{
  assert( m_ndim >= 1 );

  assert( a < m_shape[0] );

  return m_data.begin() + a*m_strides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b)
{
  assert( m_ndim >= 2 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b) const
{
  assert( m_ndim >= 2 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c)
{
  assert( m_ndim >= 3 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c) const
{
  assert( m_ndim >= 3 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d)
{
  assert( m_ndim >= 4 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d) const
{
  assert( m_ndim >= 4 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e)
{
  assert( m_ndim >= 5 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );

  return m_data.begin()+a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e) const
{
  assert( m_ndim >= 5 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );

  return m_data.begin()+a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f)
{
  assert( m_ndim >= 6 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );
  assert( f < m_shape[5] );

  return m_data.begin() +
    a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto array<X>::item(size_t a, size_t b, size_t c, size_t d, size_t e, size_t f) const
{
  assert( m_ndim >= 6 );

  assert( a < m_shape[0] );
  assert( b < m_shape[1] );
  assert( c < m_shape[2] );
  assert( d < m_shape[3] );
  assert( e < m_shape[4] );
  assert( f < m_shape[5] );

  return m_data.begin() +
    a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5];
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void array<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void array<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void array<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline array<X>& array<X>::operator*= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator/= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator+= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator-= (const array<X> &B)
{
  assert( shape() == B.shape() );
  assert( ndim()  == B.ndim()  );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X>& array<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

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
  assert( axis < m_ndim );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(this->shape(),axis), this->maxCoeff());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = this->strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), m_size);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < m_size ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::min(out[ni], m_data[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::minCoeff(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(m_ndim) );
  assert( axis >= -1 * static_cast<int>(m_ndim) );

  // get number of dimensions as integer
  int n = static_cast<int>(m_ndim);

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
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(m_ndim), true);

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
  assert( axis < m_ndim );

  // initialize output to the same shape as the input, with one axis removed
  array<X> out = array<X>::Constant(del(this->shape(),axis), this->minCoeff());

  // extended strides
  // - copy strides
  std::vector<size_t> estrides = this->strides();
  // - insert total size at the beginning
  estrides.insert(estrides.begin(), m_size);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < m_size ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] = std::max(out[ni], m_data[i]);
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::maxCoeff(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(m_ndim) );
  assert( axis >= -1 * static_cast<int>(m_ndim) );

  // get number of dimensions as integer
  int n = static_cast<int>(m_ndim);

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
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(m_ndim), true);

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

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

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
  estrides.insert(estrides.begin(), m_size);

  // extract sizes
  size_t n = estrides[axis  ];
  size_t m = estrides[axis+1];

  // perform reduction
  for ( size_t i = 0 ; i < m_size ; ++i )
  {
    // - get the new index
    size_t ni = i/n*m + i%m;
    // - store
    out[ni] += m_data[i];
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline array<X> array<X>::sum(int axis) const
{
  // check axis: (0,1,...,ndim-1) or (-1,-2,...,-ndim)
  assert( axis  <      static_cast<int>(m_ndim) );
  assert( axis >= -1 * static_cast<int>(m_ndim) );

  // get number of dimensions as integer
  int n = static_cast<int>(m_ndim);

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
  std::vector<int> axes = Private::sort_axes(axes_in, static_cast<int>(m_ndim), true);

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
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
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

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

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
inline void array<X>::abs()
{
  for ( auto &i : m_data )
    i = std::abs(i);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void array<X>::printf(std::string fmt) const
{
  if ( m_ndim == 1 )
  {
    for ( size_t j = 0 ; j < m_shape[0] ; ++j ) {
      if ( j != m_shape[0]-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
      else                     std::printf((fmt + ";\n").c_str(), (*this)(j));
    }

    return;
  }

  if ( m_ndim == 2 )
  {
    for ( size_t i = 0 ; i < m_shape[0] ; ++i ) {
      for ( size_t j = 0 ; j < m_shape[1] ; ++j ) {
        if ( j != m_shape[1]-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
        else                     std::printf((fmt + ";\n").c_str(), (*this)(i,j));
      }
    }

    return;
  }

  std::cout << "cppmat::array[";

  for ( size_t i = 0 ; i < m_ndim-1 ; ++i )
    std::cout << shape(i) << ",";

  std::cout << shape(m_ndim-1) << "]\n";

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

