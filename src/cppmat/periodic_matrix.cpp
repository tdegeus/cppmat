/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_PERIODIC_MATRIX_CPP
#define CPPMAT_PERIODIC_MATRIX_CPP

// -------------------------------------------------------------------------------------------------

#include "periodic_matrix.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline matrix<X>::matrix(const std::vector<size_t> &shape)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(shape);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>::matrix(const std::vector<size_t> &shape, X D)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(shape);

  // copy input
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline matrix<X>::matrix(const std::vector<size_t> &shape, Iterator first, Iterator last)
{
  // store shape, and other size parameters, allocate "m_data"
  resize(shape);

  // check size
  assert( m_size == last - first );

  // initialize counter
  size_t i = 0;

  // copy input
  for (auto it = first; it != last; ++it)
  {
    m_data[i] = (*it); ++i;
  }
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline void matrix<X>::resize(const std::vector<size_t> &shape)
{
  assert( shape.size()  > 0       );
  assert( shape.size() <= MAX_DIM );

  for ( size_t i = 0 ; i < MAX_DIM ; ++i )
  {
    m_shape  [i] = 1;
    m_strides[i] = 1;
  }

  m_ndim = shape.size();
  m_size = 1;

  for ( size_t i = 0 ; i < m_ndim ; ++i )
  {
    m_shape[i] = shape[i];
    m_size    *= shape[i];
  }

  for ( size_t i = 0 ; i < m_ndim ; ++i )
    for ( size_t j = i+1 ; j < m_ndim ; ++j )
      m_strides[i] *= m_shape[j];

  for ( size_t i = 0 ; i < MAX_DIM ; ++i )
    m_shape_i[i] = static_cast<int>(m_shape[i]);

  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::chdim(size_t ndim)
{
  #ifndef NDEBUG
    for ( size_t i = ndim ; i < MAX_DIM ; ++i ) assert( m_shape[i] == 1 );
  #endif

  m_ndim = ndim;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::reshape(const std::vector<size_t> &shape)
{
  #ifndef NDEBUG
    size_t n = 1;

    for ( auto &i : shape ) n *= i;

    assert( n == m_size );
  #endif

  resize(shape);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X>
inline size_t matrix<X>::size() const
{
  return m_size;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::ndim() const
{
  return m_ndim;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t matrix<X>::shape(size_t i) const
{
  assert( i < MAX_DIM );

  return m_shape[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix<X>::shape() const
{
  std::vector<size_t> ret(m_ndim);

  for ( size_t i = 0 ; i < m_ndim ; ++i ) ret[i] = m_shape[i];

  return ret;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> matrix<X>::strides(bool bytes) const
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
// index operators
// =================================================================================================

template<class X>
inline X& matrix<X>::operator[](size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator[](size_t i) const
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------


template<class X>
inline X& matrix<X>::operator()(int a)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;

  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;

  return m_data[a*m_strides[0]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;

  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;

  return m_data[a*m_strides[0]+b*m_strides[1]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b, int c)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b, int c) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b, int c, int d)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b, int c, int d) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b, int c, int d, int e)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b, int c, int d, int e) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& matrix<X>::operator()(int a, int b, int c, int d, int e, int f)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;
  f = ( f < 0 ) ? f + m_shape_i[5] : ( f >= m_shape_i[5] ) ? f - m_shape_i[5] : f ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& matrix<X>::operator()(int a, int b, int c, int d, int e, int f) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;
  f = ( f < 0 ) ? f + m_shape_i[5] : ( f >= m_shape_i[5] ) ? f - m_shape_i[5] : f ;

  return m_data[a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5]];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline X& matrix<X>::at(Iterator first, Iterator last)
{
  assert( last-first <= this->ndim() );

  int    *shape   = &m_shape  [0];
  size_t *stride  = &m_strides[0];
  size_t  idx     = 0;

  for ( auto it = first ; it != last ; ++it )
  {
    // - current index
    int i = (*it);
    // - correct for periodicity
    i = ( i < 0 ) ? i + (*shape) : ( i >= (*shape) ) ? i - (*shape) : i ;
    // - update the index
    idx += i * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return m_data[idx];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline const X& matrix<X>::at(Iterator first, Iterator last) const
{
  assert( last-first <= this->ndim() );

  int    *shape   = &m_shape  [0];
  size_t *stride  = &m_strides[0];
  size_t  idx     = 0;

  for ( auto it = first ; it != last ; ++it )
  {
    // - current index
    int i = (*it);
    // - correct for periodicity
    i = ( i < 0 ) ? i + (*shape) : ( i >= (*shape) ) ? i - (*shape) : i ;
    // - update the index
    idx += i * (*stride);
    // - move iterators forward
    ++stride;
    ++shape;
  }

  return m_data[idx];
}

// =================================================================================================
// pointer to data
// =================================================================================================

template<class X>
inline X* matrix<X>::data()
{
  return m_data.data();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X* matrix<X>::data() const
{
  return m_data.data();
}

// =================================================================================================
// iterators
// =================================================================================================

template<class X>
inline auto matrix<X>::begin()
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::begin() const
{
  return m_data.begin();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::end()
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline auto matrix<X>::end() const
{
  return m_data.end();
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;

  return m_data.begin() + a*m_strides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;

  return m_data.begin() + a*m_strides[0];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b, int c)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b, int c) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b, int c, int d)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(int a, int b, int c, int d) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;

  return m_data.begin() + a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(
  int a, int b, int c, int d, int e)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;

  return m_data.begin()+a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(
  int a, int b, int c, int d, int e) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;

  return m_data.begin()+a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(
  int a, int b, int c, int d, int e, int f)
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;
  f = ( f < 0 ) ? f + m_shape_i[5] : ( f >= m_shape_i[5] ) ? f - m_shape_i[5] : f ;

  return m_data.begin() +
    a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5];
}

// -------------------------------------------------------------------------------------------------

template<class X> inline auto matrix<X>::item(
  int a, int b, int c, int d, int e, int f) const
{
  a = ( a < 0 ) ? a + m_shape_i[0] : ( a >= m_shape_i[0] ) ? a - m_shape_i[0] : a ;
  b = ( b < 0 ) ? b + m_shape_i[1] : ( b >= m_shape_i[1] ) ? b - m_shape_i[1] : b ;
  c = ( c < 0 ) ? c + m_shape_i[2] : ( c >= m_shape_i[2] ) ? c - m_shape_i[2] : c ;
  d = ( d < 0 ) ? d + m_shape_i[3] : ( d >= m_shape_i[3] ) ? d - m_shape_i[3] : d ;
  e = ( e < 0 ) ? e + m_shape_i[4] : ( e >= m_shape_i[4] ) ? e - m_shape_i[4] : e ;
  f = ( f < 0 ) ? f + m_shape_i[5] : ( f >= m_shape_i[5] ) ? f - m_shape_i[5] : f ;

  return m_data.begin() +
    a*m_strides[0]+b*m_strides[1]+c*m_strides[2]+d*m_strides[3]+e*m_strides[4]+f*m_strides[5];
}

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void matrix<X>::arange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::zeros()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void matrix<X>::ones()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline matrix<X>& matrix<X>::operator*= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator/= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator+= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator-= (const matrix<X> &B)
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X>& matrix<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const matrix<X> &A, const matrix<X> &B)
{
  assert( A.size() == B.size() );
  assert( A.ndim() == B.ndim() );

  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const matrix<X> &A, const X &B)
{
  matrix<X> C(A.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A[i] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator* (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A * B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator/ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A / B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator+ (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A + B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline matrix<X> operator- (const X &A, const matrix<X> &B)
{
  matrix<X> C(B.shape());

  for ( size_t i = 0 ; i < C.size() ; ++i )
    C[i] = A - B[i];

  return C;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X matrix<X>::min() const
{
  return *std::min_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X matrix<X>::max() const
{
  return *std::max_element(m_data.begin(),m_data.end());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X matrix<X>::sum() const
{
  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double matrix<X>::mean() const
{
  return static_cast<double>(this->sum())/static_cast<double>(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline double matrix<X>::average(const matrix<X> &weights) const
{
  assert( size() == weights.size() );
  assert( ndim() == weights.ndim() );

  X out = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    out += m_data[i] * weights[i];

  return static_cast<double>(out)/static_cast<double>(weights.sum());
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void matrix<X>::printf(std::string fmt) const
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

  std::cout << "cppmat::periodic::matrix[";

  for ( size_t i = 0 ; i < m_ndim-1 ; ++i )
    std::cout << shape(i) << ",";

  std::cout << shape(m_ndim-1) << "]\n";

}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, matrix<X>& src)
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

  out << "cppmat::periodic::matrix[";

  for ( size_t i = 0 ; i < src.ndim()-1 ; ++i )
    out << src.shape(i) << ",";

  out << src.shape(src.ndim()-1) << "]";

  return out;
}

// =================================================================================================

}} // namespace ...

#endif

