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
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Zero(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Ones(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Constant(size_t nd, X D)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::I(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Irt(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setIrt();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Is(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setIs();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Id(size_t nd)
{
  // allocate tensor
  tensor4<X> out(nd);

  // initialize
  out.setId();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::II(size_t nd)
{
  // allocate tensor
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
  // allocate tensor
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
  // allocate tensor
  tensor2<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Zero(size_t nd)
{
  // allocate tensor
  tensor2<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Ones(size_t nd)
{
  // allocate tensor
  tensor2<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Constant(size_t nd, X D)
{
  // allocate tensor
  tensor2<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::I(size_t nd)
{
  // allocate tensor
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
  // allocate tensor
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
  // allocate tensor
  tensor2s<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Zero(size_t nd)
{
  // allocate tensor
  tensor2s<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Ones(size_t nd)
{
  // allocate tensor
  tensor2s<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Constant(size_t nd, X D)
{
  // allocate tensor
  tensor2s<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::I(size_t nd)
{
  // allocate tensor
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
  // allocate tensor
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
  // allocate tensor
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
  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d(size_t nd)
{
  // change size of internal storage
  resize(nd);

  // set dummy parameter
  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Arange(size_t nd)
{
  // allocate tensor
  tensor2d<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Zero(size_t nd)
{
  // allocate tensor
  tensor2d<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Ones(size_t nd)
{
  // allocate tensor
  tensor2d<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Constant(size_t nd, X D)
{
  // allocate tensor
  tensor2d<X> out(nd);

  // initialize
  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::I(size_t nd)
{
  // allocate tensor
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
  // allocate tensor
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
  // allocate tensor
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
  // allocate tensor
  vector<X> out(nd);

  // initialize
  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Zero(size_t nd)
{
  // allocate tensor
  vector<X> out(nd);

  // initialize
  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Ones(size_t nd)
{
  // allocate tensor
  vector<X> out(nd);

  // initialize
  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Constant(size_t nd, X D)
{
  // allocate tensor
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
  // allocate tensor
  vector<X> out(nd);

  // initialize
  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// cast to different class / type
// =================================================================================================

template<>
template<>
inline tensor2s<double> tensor2<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i ; j < m_nd ; ++j )
      out[i*m_nd-(i-1)*i/2+j-i] = ( m_data[i*m_nd+j] + m_data[j*m_nd+i] ) / 2.;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i] = m_data[i*m_nd+i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2s<double>::cast<tensor2<double>>() const
{
  tensor2<double> out(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i ; j < m_nd ; ++j )
      out[i*m_nd+j] = out[j*m_nd+i] = m_data[i*m_nd-(i-1)*i/2+j-i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2d<double> tensor2s<double>::cast<tensor2d<double>>() const
{
  tensor2d<double> out(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i] = m_data[i*m_nd-(i-1)*i/2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2<double> tensor2d<double>::cast<tensor2<double>>() const
{
  tensor2<double> out = tensor2<double>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i*m_nd+i] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline tensor2s<double> tensor2d<double>::cast<tensor2s<double>>() const
{
  tensor2s<double> out = tensor2s<double>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i*m_nd-(i-1)*i/2] = m_data[i];

  return out;
}

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2s<X>::operator tensor2<X> () const
{
  tensor2<X> out(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i ; j < m_nd ; ++j )
      out[i*m_nd+j] = out[j*m_nd+i] = m_data[i*m_nd-(i-1)*i/2+j-i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2<X> () const
{
  tensor2<double> out = tensor2<double>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i*m_nd+i] = m_data[i];

  return out;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline tensor2d<X>::operator tensor2s<X> () const
{
  tensor2s<X> out(m_nd,static_cast<X>(0));

  for ( size_t i = 0 ; i < m_nd ; ++i )
    out[i*m_nd-(i-1)*i/2] = m_data[i];

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
  m_nd   = nd;
  m_size = m_nd*m_nd*m_nd*m_nd;

  // resize storage
  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::resize(size_t nd)
{
  // store size
  m_nd   = nd;
  m_size = m_nd*m_nd;

  // resize storage
  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::resize(size_t nd)
{
  // store size
  m_nd   = nd;
  m_size = (m_nd+1)*m_nd/2;

  // resize storage
  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::resize(size_t nd)
{
  // store size
  m_nd   = nd;
  m_size = m_nd;

  // resize storage
  m_data.resize(m_size);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::resize(size_t nd)
{
  // store size
  m_nd   = nd;
  m_size = m_nd;

  // resize storage
  m_data.resize(m_size);
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X> inline size_t tensor4 <X>::size() const { return m_size; }
template<class X> inline size_t tensor2 <X>::size() const { return m_size; }
template<class X> inline size_t tensor2s<X>::size() const { return m_size; }
template<class X> inline size_t tensor2d<X>::size() const { return m_size; }
template<class X> inline size_t vector  <X>::size() const { return m_size; }
template<class X> inline size_t tensor4 <X>::ndim() const { return m_nd;   }
template<class X> inline size_t tensor2 <X>::ndim() const { return m_nd;   }
template<class X> inline size_t tensor2s<X>::ndim() const { return m_nd;   }
template<class X> inline size_t tensor2d<X>::ndim() const { return m_nd;   }
template<class X> inline size_t vector  <X>::ndim() const { return m_nd;   }

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::shape() const
{
  std::vector<size_t> out(4, m_nd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::shape() const
{
  std::vector<size_t> out(2, m_nd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> out(1, m_nd);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { m_nd*m_nd*m_nd, m_nd*m_nd, m_nd, 1 };

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
  std::vector<size_t> out = { m_nd, 1 };

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
  return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  return m_data[i*m_nd*m_nd*m_nd+j*m_nd*m_nd+k*m_nd+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2<X>::operator()(size_t i, size_t j)
{
  return m_data[i*m_nd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  return m_data[i*m_nd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X& tensor2s<X>::operator()(size_t i, size_t j)
{
  if (i <= j) return m_data[ i*m_nd - (i-1)*i/2 + j - i ];
  else        return m_data[ j*m_nd - (j-1)*j/2 + i - j ];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator()(size_t i, size_t j) const
{
  if (i <= j) return m_data[ i*m_nd - (i-1)*i/2 + j - i ];
  else        return m_data[ j*m_nd - (j-1)*j/2 + i - j ];
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

template<class X> inline X*       tensor4 <X>::data()        { return m_data.data();  }
template<class X> inline const X* tensor4 <X>::data()  const { return m_data.data();  }
template<class X> inline auto     tensor4 <X>::begin()       { return m_data.begin(); }
template<class X> inline auto     tensor4 <X>::begin() const { return m_data.begin(); }
template<class X> inline auto     tensor4 <X>::end()         { return m_data.end();   }
template<class X> inline auto     tensor4 <X>::end()   const { return m_data.end();   }
template<class X> inline X*       tensor2 <X>::data()        { return m_data.data();  }
template<class X> inline const X* tensor2 <X>::data()  const { return m_data.data();  }
template<class X> inline auto     tensor2 <X>::begin()       { return m_data.begin(); }
template<class X> inline auto     tensor2 <X>::begin() const { return m_data.begin(); }
template<class X> inline auto     tensor2 <X>::end()         { return m_data.end();   }
template<class X> inline auto     tensor2 <X>::end()   const { return m_data.end();   }
template<class X> inline X*       tensor2s<X>::data()        { return m_data.data();  }
template<class X> inline const X* tensor2s<X>::data()  const { return m_data.data();  }
template<class X> inline auto     tensor2s<X>::begin()       { return m_data.begin(); }
template<class X> inline auto     tensor2s<X>::begin() const { return m_data.begin(); }
template<class X> inline auto     tensor2s<X>::end()         { return m_data.end();   }
template<class X> inline auto     tensor2s<X>::end()   const { return m_data.end();   }
template<class X> inline X*       tensor2d<X>::data()        { return m_data.data();  }
template<class X> inline const X* tensor2d<X>::data()  const { return m_data.data();  }
template<class X> inline auto     tensor2d<X>::begin()       { return m_data.begin(); }
template<class X> inline auto     tensor2d<X>::begin() const { return m_data.begin(); }
template<class X> inline auto     tensor2d<X>::end()         { return m_data.end();   }
template<class X> inline auto     tensor2d<X>::end()   const { return m_data.end();   }
template<class X> inline X*       vector  <X>::data()        { return m_data.data();  }
template<class X> inline const X* vector  <X>::data()  const { return m_data.data();  }
template<class X> inline auto     vector  <X>::begin()       { return m_data.begin(); }
template<class X> inline auto     vector  <X>::begin() const { return m_data.begin(); }
template<class X> inline auto     vector  <X>::end()         { return m_data.end();   }
template<class X> inline auto     vector  <X>::end()   const { return m_data.end();   }

// =================================================================================================
// basic initialization
// =================================================================================================

template<class X>
inline void tensor4<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor4<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          if ( i == l and j == k )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIrt()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          if ( i == k and j == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setIs()
{
  return ( tensor4<X>::I(m_nd) + tensor4<X>::Irt(m_nd) ) / static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setId()
{
  return tensor4<X>::Is(m_nd) - tensor4<X>::II(m_nd)/static_cast<X>(m_nd);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor4<X>::setII()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          if ( i == j and k == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2s<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( m_nd * m_nd == last - first );

  // check for symmetry
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i+1 ; j < m_nd ; ++j )
      assert( first[i*m_nd+j] == first[j*m_nd+i] );
  #endif

  // copy from input (ignores lower diagonal terms)
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i ; j < m_nd ; ++j )
      m_data[i*m_nd-(i-1)*i/2+j-i] = first[i*m_nd+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopy(Iterator first, Iterator last)
{
  // check size
  assert( m_size == last - first );

  // copy
  std::copy(first, last, m_data.data());
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void tensor2d<X>::setCopyDense(Iterator first, Iterator last)
{
  // avoid compiler warning
  UNUSED(last);

  // check size
  assert( m_nd * m_nd == last - first );

  // check the input to be diagonal
  #ifndef NDEBUG
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      if ( i !=j )
        assert( !first[i*m_nd+j] );
  #endif

  // copy from input (ignores off-diagonal terms)
  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] = first[i*m_nd+i];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setArange()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(i);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setZero()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setOnes()
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setConstant(X D)
{
  for ( size_t i = 0 ; i < m_size ; ++i ) m_data[i] = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<class Iterator>
inline void vector<X>::setCopy(Iterator first, Iterator last)
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
inline tensor4<X>& tensor4<X>::operator*= (const tensor4<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const tensor4<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const tensor4<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const tensor4<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>& tensor4<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

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
inline tensor2<X>& tensor2<X>::operator*= (const tensor2<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      // - extract value
      X b = B[i*m_nd-(i-1)*i/2+j-i];
      // - store symmetrically
                    m_data[i*m_nd+j] *= b;
      if ( i != j ) m_data[j*m_nd+i] *= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      // - extract value
      X b = B[i*m_nd-(i-1)*i/2+j-i];
      // - store symmetrically
                    m_data[i*m_nd+j] /= b;
      if ( i != j ) m_data[j*m_nd+i] /= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      // - extract value
      X b = B[i*m_nd-(i-1)*i/2+j-i];
      // - store symmetrically
                    m_data[i*m_nd+j] += b;
      if ( i != j ) m_data[j*m_nd+i] += b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      // - extract value
      X b = B[i*m_nd-(i-1)*i/2+j-i];
      // - store symmetrically
                    m_data[i*m_nd+j] -= b;
      if ( i != j ) m_data[j*m_nd+i] -= b;
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = 0 ; j < m_nd ; ++j ) {
      if ( i == j ) m_data[i*m_nd+i] *= B[i];
      else          m_data[i*m_nd+j]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i*m_nd+i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const tensor2d<X> &B)
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i*m_nd+i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
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
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2s<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const tensor2s<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2s<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2s<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const tensor2d<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j=i; j<m_nd; ++j ) {
      if ( i == j ) m_data[i*m_nd-(i-1)*i/2    ] *= B[i];
      else          m_data[i*m_nd-(i-1)*i/2+j-i]  = static_cast<X>(0);
    }
  }

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const tensor2d<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i*m_nd-(i-1)*i/2] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const tensor2d<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i*m_nd-(i-1)*i/2] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
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
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2d<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const tensor2d<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const tensor2d<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2 <X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd; ++i )
    m_data[i] *= B[i*m_nd+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2 <X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd; ++i )
    m_data[i] /= B[i*m_nd+i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] *= B[i*m_nd-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const tensor2s<X> &B)
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] /= B[i*m_nd-(i-1)*i/2];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    m_data[i] -= B;

  return *this;
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
inline tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
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
inline tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
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
inline tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
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
inline vector<X>& vector<X>::operator*= (const vector<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const vector<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const vector<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const vector<X> &B)
{
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B[i];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator*= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] *= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator/= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator+= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] += B;

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator-= (const X &B)
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] -= B;

  return *this;
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

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          for ( size_t m = 0 ; m < m_nd ; ++m )
            for ( size_t n = 0 ; n < m_nd ; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l) * B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
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

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
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

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += ( m_data[i*m_nd-(i-1)*i/2] * B[i*m_nd-(i-1)*i/2] );

  for ( size_t i = 0 ; i<m_nd ; ++i )
    for ( size_t j = i+1 ; j<m_nd ; ++j )
      C += ( static_cast<X>(2) * m_data[i*m_nd-(i-1)*i/2+j-i] * B[i*m_nd-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += m_data[i*m_nd-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      for ( size_t l = 0 ; l < m_nd ; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += m_data[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += m_data[i]*B[i*m_nd-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += m_data[i]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2d<X> C = tensor2d<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
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

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      for ( size_t l = 0 ; l < m_nd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      for ( size_t l = 0 ; l < m_nd ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t k = 0 ; k < m_nd ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::cross(const vector<X> &B) const
{
  assert( size() == B.size() );
  assert( ndim() == B.ndim() );

  if ( m_nd != 3 )
    throw std::runtime_error("'cross' only implemented in 3D");

  vector<X> C(3);

  C[0] =                     m_data[1]*B[2]-B[1]*m_data[2] ;
  C[1] = static_cast<X>(-1)*(m_data[0]*B[2]-B[0]*m_data[2]);
  C[2] =                     m_data[0]*B[1]-B[0]*m_data[1] ;

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
  tensor4<X> C(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::RT() const
{
  tensor4<X> C(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::LT() const
{
  tensor4<X> C(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::T() const
{
  tensor2<X> C(m_nd);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      C(j,i) = (*this)(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::T() const
{
  tensor2s<X> C(m_nd);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C[i] = (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::T() const
{
  tensor2d<X> C(m_nd);

  for ( size_t i = 0 ; i < m_size ; ++i )
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

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_nd ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::det() const
{
  if ( m_nd==2 )
   return m_data[0] * m_data[3] - m_data[1] * m_data[2];

  if ( m_nd==3 )
    return ( m_data[0] * m_data[4] * m_data[8] +
             m_data[1] * m_data[5] * m_data[6] +
             m_data[2] * m_data[3] * m_data[7] ) -
           ( m_data[2] * m_data[4] * m_data[6] +
             m_data[1] * m_data[3] * m_data[8] +
             m_data[0] * m_data[5] * m_data[7] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::det() const
{
  if ( m_nd==2 )
   return m_data[0] * m_data[2] - m_data[1] * m_data[1];

  if ( m_nd==3 )
    return (                     m_data[0] * m_data[3] * m_data[5] +
             static_cast<X>(2) * m_data[1] * m_data[2] * m_data[4] ) -
           (                     m_data[4] * m_data[4] * m_data[0] +
                                 m_data[2] * m_data[2] * m_data[3] +
                                 m_data[1] * m_data[1] * m_data[5] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::det() const
{
  X C = static_cast<X>(1);

  for ( size_t i = 0 ; i < m_nd ; ++i )
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
  tensor2<X> C(m_nd);

  if ( m_nd==2 )
  {
    C[0] =                      m_data[3] / D;
    C[1] = static_cast<X>(-1) * m_data[1] / D;
    C[2] = static_cast<X>(-1) * m_data[2] / D;
    C[3] =                      m_data[0] / D;

    return C;
  }

  if ( m_nd==3 )
  {
    C[0] = (m_data[4]*m_data[8]-m_data[5]*m_data[7]) / D;
    C[1] = (m_data[2]*m_data[7]-m_data[1]*m_data[8]) / D;
    C[2] = (m_data[1]*m_data[5]-m_data[2]*m_data[4]) / D;
    C[3] = (m_data[5]*m_data[6]-m_data[3]*m_data[8]) / D;
    C[4] = (m_data[0]*m_data[8]-m_data[2]*m_data[6]) / D;
    C[5] = (m_data[2]*m_data[3]-m_data[0]*m_data[5]) / D;
    C[6] = (m_data[3]*m_data[7]-m_data[4]*m_data[6]) / D;
    C[7] = (m_data[1]*m_data[6]-m_data[0]*m_data[7]) / D;
    C[8] = (m_data[0]*m_data[4]-m_data[1]*m_data[3]) / D;

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
  tensor2s<X> C(m_nd);

  if ( m_nd==2 )
  {
    C[0] =                      m_data[2] / D;
    C[1] = static_cast<X>(-1) * m_data[1] / D;
    C[2] =                      m_data[0] / D;

    return C;
  }

  if ( m_nd==3 )
  {
    C[0] = (m_data[3]*m_data[5]-m_data[4]*m_data[4]) / D;
    C[1] = (m_data[2]*m_data[4]-m_data[1]*m_data[5]) / D;
    C[2] = (m_data[1]*m_data[4]-m_data[2]*m_data[3]) / D;
    C[3] = (m_data[0]*m_data[5]-m_data[2]*m_data[2]) / D;
    C[4] = (m_data[2]*m_data[1]-m_data[0]*m_data[4]) / D;
    C[5] = (m_data[0]*m_data[3]-m_data[1]*m_data[1]) / D;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C(m_nd);

  for ( size_t i = 0; i < m_nd ; ++i )
    C[i] = static_cast<X>(1) / m_data[i];

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
  assert( m_size == B.size() );
  assert( m_nd   == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2s<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      if ( m_data[i*m_nd+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2d<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      if ( m_data[i*m_nd+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(i,j) ) return false;
      if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(j,i) ) return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i ; j < m_nd ; ++j )
      if ( m_data[i*m_nd-(i-1)*i/2+j-i] != B(i,j) ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2d<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_size ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = 0 ; j < m_nd ; ++j ) {
      if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
      else          { if (              B(i,j) ) return false; }
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2s<X> &B) const
{
  assert( m_nd == B.ndim() );

  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = i ; j < m_nd ; ++j ) {
      if ( i == j ) { if ( m_data[i] != B(i,i) ) return false; }
      else          { if (              B(i,j) ) return false; }
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const vector<X> &B) const
{
  for ( size_t i = 0 ; i < m_size ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// =================================================================================================
// structure check
// =================================================================================================

template<class X>
inline bool tensor2<X>::issymmetric() const
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i+1 ; j < m_nd ; ++j )
      if ( m_data[ i*m_nd + j ] != m_data[ j*m_nd + i ] )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      if ( i != j )
        if ( m_data[ i*m_nd + j ] )
          return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = i+1 ; j < m_nd ; ++j )
      if ( m_data[i*m_nd-(i-1)*i/2+j-i] )
        return false;

  return true;
}

// =================================================================================================
// basic algebra
// =================================================================================================

template<class X>
inline X tensor4<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::abs(m_data[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::abs(m_data[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::abs(m_data[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::abs(m_data[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::abs(m_data[i]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::length() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < m_size ; ++i )
    C += std::pow(m_data[i],2.);

  return std::sqrt(C);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::setUnitLength()
{
  X C = length();

  if ( C <= static_cast<X>(0) ) return;

  for ( size_t i = 0 ; i < m_size ; ++i )
    m_data[i] /= C;
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void tensor4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(m_nd).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < m_nd ; ++i )
    for ( size_t j = 0 ; j < m_nd ; ++j )
      for ( size_t k = 0 ; k < m_nd ; ++k )
        for ( size_t l = 0 ; l < m_nd ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, (*this)(i,j,k,l) );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = 0 ; j < m_nd ; ++j ) {
      if ( j != m_nd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else               std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = 0 ; j < m_nd ; ++j ) {
      if ( j != m_nd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else               std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::printf(std::string fmt) const
{
  for ( size_t i = 0 ; i < m_nd ; ++i ) {
    for ( size_t j = 0 ; j < m_nd ; ++j ) {
      if ( j != m_nd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(i,j));
      else               std::printf((fmt + ";\n").c_str(), (*this)(i,j));
    }
  }
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  for ( size_t j = 0 ; j < m_nd ; ++j ) {
    if ( j != m_nd-1 ) std::printf((fmt + ","  ).c_str(), (*this)(j));
    else               std::printf((fmt + ";\n").c_str(), (*this)(j));
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

