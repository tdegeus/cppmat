/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR2_CPP
#define CPPMAT_TENSOR2_CPP

// -------------------------------------------------------------------------------------------------

#include "tensor2.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian2d {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::Arange()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Zero()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Ones()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Constant(X D)
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::I()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Irt()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setIrt();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Is()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setIs();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Id()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setId();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::II()
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setII();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor4<X> tensor4<X>::Copy(Iterator first)
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setCopy(first,first+16);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor4<X> tensor4<X>::Copy(Iterator first, Iterator last)
{
  // allocate tensor
  tensor4<X> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor4<X>::Size()
{
  return 16;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Arange()
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Zero()
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Ones()
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Constant(X D)
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::I()
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2<X> tensor2<X>::Copy(Iterator first)
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setCopy(first,first+4);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2<X> tensor2<X>::Copy(Iterator first, Iterator last)
{
  // allocate tensor
  tensor2<X> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2<X>::Size()
{
  return 4;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Arange()
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Zero()
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Ones()
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Constant(X D)
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::I()
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2s<X> tensor2s<X>::Copy(Iterator first)
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setCopy(first,first+3);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2s<X> tensor2s<X>::Copy(Iterator first, Iterator last)
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2s<X> tensor2s<X>::CopyDense(Iterator first, Iterator last)
{
  // allocate tensor
  tensor2s<X> out;

  // initialize
  out.setCopyDense(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2s<X>::Size()
{
  return 3;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d()
{
  // set dummy parameter
  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Arange()
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Zero()
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Ones()
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Constant(X D)
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::I()
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2d<X> tensor2d<X>::Copy(Iterator first)
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setCopy(first,first+2);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2d<X> tensor2d<X>::Copy(Iterator first, Iterator last)
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline tensor2d<X> tensor2d<X>::CopyDense(Iterator first, Iterator last)
{
  // allocate tensor
  tensor2d<X> out;

  // initialize
  out.setCopyDense(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t tensor2d<X>::Size()
{
  return 2;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Arange()
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Zero()
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Ones()
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Constant(X D)
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline vector<X> vector<X>::Copy(Iterator first)
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setCopy(first,first+2);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline vector<X> vector<X>::Copy(Iterator first, Iterator last)
{
  // allocate tensor
  vector<X> out;

  // initialize
  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline size_t vector<X>::Size()
{
  return 2;
}


// =================================================================================================
// cast to different class / type
// =================================================================================================

template<>
template<>
inline tensor2s<double> tensor2<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out;

  out[0] =   m_data[0];
  out[1] = ( m_data[1] + m_data[2] ) / 2.;
  out[2] =   m_data[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2s<double>::cast<tensor2<double>>() const
{
  tensor2<double> out;

  out[0]          = m_data[0];
  out[1] = out[2] = m_data[1];
  out[3]          = m_data[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2s<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2d<double>::cast<tensor2<double>>() const
{
  tensor2<double> out;

  out[0] = m_data[0];
  out[3] = m_data[1];

  out[1] = out[2] = 0.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2s<double> tensor2d<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out;

  out[0] = m_data[0];
  out[2] = m_data[1];

  out[1] = 0.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2s<X>::operator tensor2<X> () const
{
  tensor2<X> out;

  out[0]          = m_data[0];
  out[1] = out[2] = m_data[1];
  out[3]          = m_data[2];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2<X> () const
{
  tensor2<X> out;

  out[0] = m_data[0];
  out[3] = m_data[1];

  out[1] = out[2] = static_cast<X>(0);

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2s<X> () const
{
  tensor2s<X> out;

  out[0] = m_data[0];
  out[2] = m_data[1];

  out[1] = static_cast<X>(0);

  return out;
}
#endif

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X> inline size_t tensor4 <X>::size() const { return 16; }
template<class X> inline size_t tensor2 <X>::size() const { return 4;  }
template<class X> inline size_t tensor2s<X>::size() const { return 3;  }
template<class X> inline size_t tensor2d<X>::size() const { return 2;  }
template<class X> inline size_t vector  <X>::size() const { return 2;  }
template<class X> inline size_t tensor4 <X>::ndim() const { return 2;  }
template<class X> inline size_t tensor2 <X>::ndim() const { return 2;  }
template<class X> inline size_t tensor2s<X>::ndim() const { return 2;  }
template<class X> inline size_t tensor2d<X>::ndim() const { return 2;  }
template<class X> inline size_t vector  <X>::ndim() const { return 2;  }

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::shape() const
{
  std::vector<size_t> out(4,2);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::shape() const
{
  std::vector<size_t> out(2,2);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> out(1,2);

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
// index operators
// =================================================================================================

template<class X> inline X&       tensor4 <X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& tensor4 <X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline X&       tensor2 <X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& tensor2 <X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline X&       tensor2s<X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& tensor2s<X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline X&       tensor2d<X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& tensor2d<X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline X&       vector  <X>::operator[](size_t i)       { return m_data[i]; }
template<class X> inline const X& vector  <X>::operator[](size_t i) const { return m_data[i]; }

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l)
{
  return m_data[i*8+j*4+k*2+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  return m_data[i*8+j*4+k*2+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2<X>::operator()(size_t i, size_t j)
{
  return m_data[i*2+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  return m_data[i*2+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2s<X>::operator()(size_t i, size_t j)
{
  if ( i == 0 ) {
    if ( j == 0 ) return m_data[0];
    else          return m_data[1];
  }
  else {
    if ( j == 0 ) return m_data[1];
    else          return m_data[2];
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator()(size_t i, size_t j) const
{
  if ( i == 0 ) {
    if ( j == 0 ) return m_data[0];
    else          return m_data[1];
  }
  else {
    if ( j == 0 ) return m_data[1];
    else          return m_data[2];
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2d<X>::operator()(size_t i, size_t j)
{
  if (i == j) return m_data[i];
  else        return m_zero[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2d<X>::operator()(size_t i, size_t j) const
{
  if (i == j) return m_data[i];
  else        return m_zero[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& vector<X>::operator()(size_t i)
{
  return m_data[i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& vector<X>::operator()(size_t i) const
{
  return m_data[i];
}

// =================================================================================================
// pointers / iterators
// =================================================================================================

template<class X> inline X*       tensor4 <X>::data()        { return std::begin(m_data);      }
template<class X> inline const X* tensor4 <X>::data()  const { return std::begin(m_data);      }
template<class X> inline auto     tensor4 <X>::begin()       { return std::begin(m_data);      }
template<class X> inline auto     tensor4 <X>::begin() const { return std::begin(m_data);      }
template<class X> inline auto     tensor4 <X>::end()         { return std::begin(m_data) + 16; }
template<class X> inline auto     tensor4 <X>::end()   const { return std::begin(m_data) + 16; }
template<class X> inline X*       tensor2 <X>::data()        { return std::begin(m_data);      }
template<class X> inline const X* tensor2 <X>::data()  const { return std::begin(m_data);      }
template<class X> inline auto     tensor2 <X>::begin()       { return std::begin(m_data);      }
template<class X> inline auto     tensor2 <X>::begin() const { return std::begin(m_data);      }
template<class X> inline auto     tensor2 <X>::end()         { return std::begin(m_data) + 4;  }
template<class X> inline auto     tensor2 <X>::end()   const { return std::begin(m_data) + 4;  }
template<class X> inline X*       tensor2s<X>::data()        { return std::begin(m_data);      }
template<class X> inline const X* tensor2s<X>::data()  const { return std::begin(m_data);      }
template<class X> inline auto     tensor2s<X>::begin()       { return std::begin(m_data);      }
template<class X> inline auto     tensor2s<X>::begin() const { return std::begin(m_data);      }
template<class X> inline auto     tensor2s<X>::end()         { return std::begin(m_data) + 3;  }
template<class X> inline auto     tensor2s<X>::end()   const { return std::begin(m_data) + 3;  }
template<class X> inline X*       tensor2d<X>::data()        { return std::begin(m_data);      }
template<class X> inline const X* tensor2d<X>::data()  const { return std::begin(m_data);      }
template<class X> inline auto     tensor2d<X>::begin()       { return std::begin(m_data);      }
template<class X> inline auto     tensor2d<X>::begin() const { return std::begin(m_data);      }
template<class X> inline auto     tensor2d<X>::end()         { return std::begin(m_data) + 2;  }
template<class X> inline auto     tensor2d<X>::end()   const { return std::begin(m_data) + 2;  }
template<class X> inline X*       vector  <X>::data()        { return std::begin(m_data);      }
template<class X> inline const X* vector  <X>::data()  const { return std::begin(m_data);      }
template<class X> inline auto     vector  <X>::begin()       { return std::begin(m_data);      }
template<class X> inline auto     vector  <X>::begin() const { return std::begin(m_data);      }
template<class X> inline auto     vector  <X>::end()         { return std::begin(m_data) + 2;  }
template<class X> inline auto     vector  <X>::end()   const { return std::begin(m_data) + 2;  }

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void tensor4<X>::setArange()
{
  for ( size_t i = 0 ; i < 16 ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setZero()
{
  for ( size_t i = 0 ; i < 16 ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setOnes()
{
  for ( size_t i = 0 ; i < 16 ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < 16 ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor4<X>::setCopy(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 16 == last - first );

  // copy input
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] = first[i  ];
    m_data[i+1] = first[i+1];
    m_data[i+2] = first[i+2];
    m_data[i+3] = first[i+3];
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          if ( i == l and j == k )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIrt()
{
  this->setZero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          if ( i == k and j == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIs()
{
  return ( tensor4<X>::I() + tensor4<X>::Irt() ) / static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setId()
{
  return tensor4<X>::Is() - tensor4<X>::II()/static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setII()
{
  this->setZero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          if ( i == j and k == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setArange()
{
  for ( size_t i = 0 ; i < 4 ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setZero()
{
  m_data[0] = m_data[1] = m_data[2] = m_data[3] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setOnes()
{
  m_data[0] = m_data[1] = m_data[2] = m_data[3] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setConstant(X D)
{
  m_data[0] = m_data[1] = m_data[2] = m_data[3] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2<X>::setCopy(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 4 == last - first );

  // copy input
  m_data[0] = first[0];
  m_data[1] = first[1];
  m_data[2] = first[2];
  m_data[3] = first[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setI()
{
  m_data[0] = static_cast<X>(1);
  m_data[1] = static_cast<X>(0);
  m_data[2] = static_cast<X>(0);
  m_data[3] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setArange()
{
  for ( size_t i = 0 ; i < 3 ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setZero()
{
  m_data[0] = m_data[1] = m_data[2] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setOnes()
{
  m_data[0] = m_data[1] = m_data[2] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setConstant(X D)
{
  m_data[0] = m_data[1] = m_data[2] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopy(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 3 == last - first );

  // copy input
  m_data[0] = first[0];
  m_data[1] = first[1];
  m_data[2] = first[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 4 == last - first );

  // check for symmetry
  assert( first[1] == first[2] );

  // copy from input (ignores lower diagonal terms)
  m_data[0] = first[0];
  m_data[1] = first[1];
  m_data[2] = first[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setI()
{
  m_data[0] = static_cast<X>(1);
  m_data[1] = static_cast<X>(0);
  m_data[2] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setArange()
{
  for ( size_t i = 0 ; i < 2 ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setZero()
{
  m_data[0] = m_data[1] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setOnes()
{
  m_data[0] = m_data[1] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setConstant(X D)
{
  m_data[0] = m_data[1] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopy(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 2 == last - first );

  // copy input
  m_data[0] = first[0];
  m_data[1] = first[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 4 == last - first );

  // check the input to be diagonal
  assert( ! first[1] );
  assert( ! first[2] );

  // copy from input (ignores off-diagonal terms)
  m_data[0] = first[0];
  m_data[1] = first[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setI()
{
  m_data[0] = m_data[1] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setArange()
{
  for ( size_t i = 0 ; i < 2 ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setZero()
{
  m_data[0] = m_data[1] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setOnes()
{
  m_data[0] = m_data[1] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setConstant(X D)
{
  m_data[0] = m_data[1] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void vector<X>::setCopy(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( 2 == last - first );

  // copy input
  m_data[0] = first[0];
  m_data[1] = first[1];
}

// =================================================================================================
// arithmetic operators - tensor4
// =================================================================================================

template<class X>
inline tensor4<X>& tensor4<X>::operator*= (const tensor4<X> &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] *= B[i  ];
    m_data[i+1] *= B[i+1];
    m_data[i+2] *= B[i+2];
    m_data[i+3] *= B[i+3];
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const tensor4<X> &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] /= B[i  ];
    m_data[i+1] /= B[i+1];
    m_data[i+2] /= B[i+2];
    m_data[i+3] /= B[i+3];
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const tensor4<X> &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] += B[i  ];
    m_data[i+1] += B[i+1];
    m_data[i+2] += B[i+2];
    m_data[i+3] += B[i+3];
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const tensor4<X> &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] -= B[i  ];
    m_data[i+1] -= B[i+1];
    m_data[i+2] -= B[i+2];
    m_data[i+3] -= B[i+3];
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] *= B;
    m_data[i+1] *= B;
    m_data[i+2] *= B;
    m_data[i+3] *= B;
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] /= B;
    m_data[i+1] /= B;
    m_data[i+2] /= B;
    m_data[i+3] /= B;
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] += B;
    m_data[i+1] += B;
    m_data[i+2] += B;
    m_data[i+3] += B;
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    m_data[i  ] -= B;
    m_data[i+1] -= B;
    m_data[i+2] -= B;
    m_data[i+3] -= B;
  }

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2
// =================================================================================================

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[1];
  m_data[2] *= B[2];
  m_data[3] *= B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2<X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[1];
  m_data[2] /= B[2];
  m_data[3] /= B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2<X> &B)
{
  m_data[0] += B[0];
  m_data[1] += B[1];
  m_data[2] += B[2];
  m_data[3] += B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2<X> &B)
{
  m_data[0] -= B[0];
  m_data[1] -= B[1];
  m_data[2] -= B[2];
  m_data[3] -= B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2s<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[1]; m_data[2] *= B[1];
  m_data[3] *= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2s<X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[1]; m_data[2] /= B[1];
  m_data[3] /= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2s<X> &B)
{
  m_data[0] += B[0];
  m_data[1] += B[1]; m_data[2] += B[1];
  m_data[3] += B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2s<X> &B)
{
  m_data[0] -= B[0];
  m_data[1] -= B[1]; m_data[2] -= B[1];
  m_data[3] -= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2d<X> &B)
{
  m_data[0] *= B[0];
  m_data[3] *= B[1];
  m_data[1] = m_data[2] = static_cast<X>(0);

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2d<X> &B)
{
  m_data[0] += B[0];
  m_data[3] += B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2d<X> &B)
{
  m_data[0] -= B[0];
  m_data[3] -= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const X &B)
{
  m_data[0] *= B;
  m_data[1] *= B;
  m_data[2] *= B;
  m_data[3] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const X &B)
{
  m_data[0] /= B;
  m_data[1] /= B;
  m_data[2] /= B;
  m_data[3] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const X &B)
{
  m_data[0] += B;
  m_data[1] += B;
  m_data[2] += B;
  m_data[3] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const X &B)
{
  m_data[0] -= B;
  m_data[1] -= B;
  m_data[2] -= B;
  m_data[3] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2s
// =================================================================================================

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2s<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[1];
  m_data[2] *= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const tensor2s<X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[1];
  m_data[2] /= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2s<X> &B)
{
  m_data[0] += B[0];
  m_data[1] += B[1];
  m_data[2] += B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2s<X> &B)
{
  m_data[0] -= B[0];
  m_data[1] -= B[1];
  m_data[2] -= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2d<X> &B)
{
  m_data[0] *= B[0];
  m_data[2] *= B[1];
  m_data[1]  = static_cast<X>(0);

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2d<X> &B)
{
  m_data[0] += B[0];
  m_data[2] += B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2d<X> &B)
{
  m_data[0] -= B[0];
  m_data[2] -= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const X &B)
{
  m_data[0] *= B;
  m_data[1] *= B;
  m_data[2] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const X &B)
{
  m_data[0] /= B;
  m_data[1] /= B;
  m_data[2] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const X &B)
{
  m_data[0] += B;
  m_data[1] += B;
  m_data[2] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const X &B)
{
  m_data[0] -= B;
  m_data[1] -= B;
  m_data[2] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - tensor2d
// =================================================================================================

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2d<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const tensor2d<X> &B)
{
  m_data[0] += B[0];
  m_data[1] += B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const tensor2d<X> &B)
{
  m_data[0] -= B[0];
  m_data[1] -= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2 <X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2 <X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[3];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2s<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2s<X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const X &B)
{
  m_data[0] *= B;
  m_data[1] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const X &B)
{
  m_data[0] /= B;
  m_data[1] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const X &B)
{
  m_data[0] += B;
  m_data[1] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const X &B)
{
  m_data[0] -= B;
  m_data[1] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - vector
// =================================================================================================

template<class X>
inline vector<X>& vector<X>::operator*= (const vector<X> &B)
{
  m_data[0] *= B[0];
  m_data[1] *= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  m_data[0] /= B[0];
  m_data[1] /= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  m_data[0] += B[0];
  m_data[1] += B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const vector<X> &B)
{
  m_data[0] -= B[0];
  m_data[1] -= B[1];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator*= (const X &B)
{
  m_data[0] *= B;
  m_data[1] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const X &B)
{
  m_data[0] /= B;
  m_data[1] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const X &B)
{
  m_data[0] += B;
  m_data[1] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const X &B)
{
  m_data[0] -= B;
  m_data[1] -= B;

  return *this;
}

// =================================================================================================
// arithmetic operators - mixed class
// =================================================================================================

template<class X>
inline tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[2] * B[1];
  C[3] = A[3] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[2] / B[1];
  C[3] = A[3] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[2] + B[1];
  C[3] = A[3] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[2] - B[1];
  C[3] = A[3] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[1] * B[2];
  C[3] = A[2] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[1] / B[2];
  C[3] = A[2] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[1] + B[2];
  C[3] = A[2] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[1] - B[2];
  C[3] = A[2] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] = A[1] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] = A[1] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] = A[1] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[2] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const X &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];

  return C;
}

// =================================================================================================
// arithmetic operators - mixed class - view ? normal
// =================================================================================================

template<class X>
inline tensor4<X> operator* (const map::tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const map::tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const map::tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const map::tensor4<X> &A, const tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator* (const map::tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const map::tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const map::tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const map::tensor4<X> &A, const X &B)
{
  tensor4<X> C;

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
inline tensor2<X> operator* (const map::tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const map::tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const map::tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[2] * B[1];
  C[3] = A[3] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const map::tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[2] / B[1];
  C[3] = A[3] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[2] + B[1];
  C[3] = A[3] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2<X> &A, const tensor2s<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[2] - B[1];
  C[3] = A[3] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2<X> &A, const tensor2d<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const map::tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const map::tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2<X> &A, const X &B)
{
  tensor2<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const map::tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[1] * B[2];
  C[3] = A[2] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const map::tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[1] / B[2];
  C[3] = A[2] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[1] + B[2];
  C[3] = A[2] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2s<X> &A, const tensor2<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[1] - B[2];
  C[3] = A[2] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const map::tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const map::tensor2d<X> &A, const tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const map::tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const map::tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const map::tensor2s<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const map::tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const map::tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const map::tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const map::tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const map::tensor2s<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const map::tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const map::tensor2d<X> &A, const X &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] = A[1] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const map::tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] = A[1] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const map::tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] = A[1] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator+ (const map::tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator- (const map::tensor2d<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const map::tensor2d<X> &A, const X &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2 <X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const map::tensor2s<X> &A, const tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[2] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const map::vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const map::vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const map::vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const map::vector<X> &A, const vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const map::vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const map::vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const map::vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const map::vector<X> &A, const X &B)
{
  vector<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;

  return C;
}















// =================================================================================================
// arithmetic operators - mixed class - normal ? view
// =================================================================================================

template<class X>
inline tensor4<X> operator* (const tensor4<X> &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const tensor4<X> &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const tensor4<X> &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const tensor4<X> &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator* (const X &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator/ (const X &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator+ (const X &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor4<X> operator- (const X &A, const map::tensor4<X> &B)
{
  tensor4<X> C;

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
inline tensor2<X> operator* (const tensor2<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2<X> &A, const map::tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[2] * B[1];
  C[3] = A[3] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2<X> &A, const map::tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[2] / B[1];
  C[3] = A[3] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const map::tensor2s<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[2] + B[1];
  C[3] = A[3] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const map::tensor2s<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[2] - B[1];
  C[3] = A[3] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2<X> &A, const map::tensor2d<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2<X> &A, const map::tensor2d<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const tensor2s<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[2] = A[1] * B[2];
  C[3] = A[2] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const tensor2s<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[2] = A[1] / B[2];
  C[3] = A[2] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2s<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[2] = A[1] + B[2];
  C[3] = A[2] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2s<X> &A, const map::tensor2<X> &B)
{
  tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[2] = A[1] - B[2];
  C[3] = A[2] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const tensor2d<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const tensor2d<X> &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator* (const X &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator/ (const X &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator+ (const X &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> operator- (const X &A, const map::tensor2<X> &B)
{
  tensor2<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const tensor2s<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const tensor2s<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2s<X> &A, const map::tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2s<X> &A, const map::tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const tensor2d<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] = A[1] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const tensor2d<X> &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] = A[1] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator* (const X &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator/ (const X &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const map::tensor2s<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator+ (const X &A, const map::tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> operator- (const X &A, const map::tensor2d<X> &B)
{
  tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator+ (const tensor2d<X> &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator- (const tensor2d<X> &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const map::tensor2 <X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2d<X> &A, const map::tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator/ (const tensor2d<X> &A, const map::tensor2s<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2 <X> &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const tensor2s<X> &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[2] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> operator* (const X &A, const map::tensor2d<X> &B)
{
  tensor2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const vector<X> &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const vector<X> &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const vector<X> &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const vector<X> &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator* (const X &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator/ (const X &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator+ (const X &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> operator- (const X &A, const map::vector<X> &B)
{
  vector<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];

  return C;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          for ( size_t m = 0 ; m < 2 ; ++m )
            for ( size_t n = 0 ; n < 2 ; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l) * B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[3]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  X C;

  C  = m_data[0] * B[0];
  C += m_data[1] * B[1] * static_cast<X>(2);
  C += m_data[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[2]*B[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      for ( size_t l = 0 ; l < 2 ; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[1];
}


// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  tensor2d<X> C = tensor2d<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  vector<X> C = vector<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    C(i) += (*this)(i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::dot(const vector<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 2 ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      for ( size_t l = 0 ; l < 2 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      for ( size_t l = 0 ; l < 2 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  tensor4<X> C = tensor4<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  tensor2<X> C = tensor2<X>::Zero();

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i) * B(j);

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

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline tensor4<X> tensor4<X>::T() const
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::RT() const
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::LT() const
{
  tensor4<X> C;

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::T() const
{
  tensor2<X> C;

  C[0] = m_data[0];
  C[2] = m_data[1];
  C[1] = m_data[2];
  C[3] = m_data[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::T() const
{
  tensor2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::T() const
{
  tensor2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];

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
  return m_data[0] + m_data[3];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::trace() const
{
  return m_data[0] + m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::trace() const
{
  return m_data[0] + m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::det() const
{
  return m_data[0] * m_data[3] - m_data[1] * m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::det() const
{
  return m_data[0] * m_data[2] - m_data[1] * m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::det() const
{
  return m_data[0] * m_data[1];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2<X> C;

  C[0] =                      m_data[3] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] = static_cast<X>(-1) * m_data[2] / D;
  C[3] =                      m_data[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2s<X> C;

  C[0] =                      m_data[2] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] =                      m_data[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C;

  C[0] = static_cast<X>(1) / m_data[0];
  C[1] = static_cast<X>(1) / m_data[1];

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
  for ( size_t i = 0 ; i < 16 ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;
  if ( m_data[3] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[1] ) return false;
  if ( m_data[3] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2d<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[3] != B[1] ) return false;
  if ( m_data[1]         ) return false;
  if ( m_data[2]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[1] != B[2] ) return false;
  if ( m_data[2] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2d<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[2] != B[1] ) return false;
  if ( m_data[1]         ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2d<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[3] ) return false;
  if (              B[1] ) return false;
  if (              B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[2] ) return false;
  if (              B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const vector<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;

  return true;
}

// =================================================================================================
// structure check
// =================================================================================================

template<class X>
inline bool tensor2<X>::issymmetric() const
{
  if ( m_data[1] != m_data[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::isdiagonal() const
{
  if ( m_data[1] != static_cast<X>(0) ) return false;
  if ( m_data[2] != static_cast<X>(0) ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::isdiagonal() const
{
  if ( m_data[1] != static_cast<X>(0) ) return false;

  return true;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X tensor4<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 16 ; i += 4 )
  {
    C += std::abs(m_data[i  ]);
    C += std::abs(m_data[i+1]);
    C += std::abs(m_data[i+2]);
    C += std::abs(m_data[i+3]);
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);
  C += std::abs(m_data[2]);
  C += std::abs(m_data[3]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);
  C += std::abs(m_data[2]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::length() const
{
  X C;

  C  = std::pow(m_data[0],2.);
  C += std::pow(m_data[1],2.);

  return std::sqrt(C);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setUnitLength()
{
  X C = length();

  if ( C <= static_cast<X>(0) ) return;

  m_data[0] /= C;
  m_data[1] /= C;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void tensor4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(2).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
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
  std::printf((fmt + "," + fmt + "\n").c_str(), m_data[0], m_data[1]);
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

}} // namespace ...

#endif

