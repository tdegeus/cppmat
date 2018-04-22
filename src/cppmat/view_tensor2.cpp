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
inline tensor4<X>::tensor4()
{
  m_data = nullptr;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X> tensor4<X>::Map(const X *D)
{
  // allocate matrix
  tensor4<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>::tensor2()
{
  m_data = nullptr;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X> tensor2<X>::Map(const X *D)
{
  // allocate matrix
  tensor2<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>::tensor2s()
{
  m_data = nullptr;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X> tensor2s<X>::Map(const X *D)
{
  // allocate matrix
  tensor2s<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d()
{
  m_data = nullptr;

  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X> tensor2d<X>::Map(const X *D)
{
  // allocate matrix
  tensor2d<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector()
{
  m_data = nullptr;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X> vector<X>::Map(const X *D)
{
  // allocate matrix
  vector<X> out;

  // initialize
  out.setMap(D);

  return out;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X> inline void tensor4 <X>::setMap(const X *D) { m_data = D; }
template<class X> inline void tensor2 <X>::setMap(const X *D) { m_data = D; }
template<class X> inline void tensor2s<X>::setMap(const X *D) { m_data = D; }
template<class X> inline void tensor2d<X>::setMap(const X *D) { m_data = D; }
template<class X> inline void vector  <X>::setMap(const X *D) { m_data = D; }

// =================================================================================================
// cast to different class / type
// =================================================================================================

template<>
template<>
inline reg::tensor2s<double> tensor2<double>::cast<reg::tensor2s<double>>() const
{
  reg::tensor2s<double> out;

  out[0] =   m_data[0];
  out[1] = ( m_data[1] + m_data[2] ) / 2.;
  out[2] =   m_data[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2d<double> tensor2<double>::cast<reg::tensor2d<double>>() const
{
  reg::tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[3];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2<double> tensor2s<double>::cast<reg::tensor2<double>>() const
{
  reg::tensor2<double> out;

  out[0]          = m_data[0];
  out[1] = out[2] = m_data[1];
  out[3]          = m_data[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2d<double> tensor2s<double>::cast<reg::tensor2d<double>>() const
{
  reg::tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[2];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2<double> tensor2d<double>::cast<reg::tensor2<double>>() const
{
  reg::tensor2<double> out;

  out[0] = m_data[0];
  out[3] = m_data[1];

  out[1] = out[2] = 0.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline reg::tensor2s<double> tensor2d<double>::cast<reg::tensor2s<double>>() const
{
  reg::tensor2s<double> out;

  out[0] = m_data[0];
  out[2] = m_data[1];

  out[1] = 0.0;

  return out;
}

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

template<class X> inline const X& tensor4 <X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline const X& tensor2 <X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline const X& tensor2s<X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline const X& tensor2d<X>::operator[](size_t i) const { return m_data[i]; }
template<class X> inline const X& vector  <X>::operator[](size_t i) const { return m_data[i]; }

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor4<X>::operator()(size_t i, size_t j, size_t k, size_t l) const
{
  return m_data[i*8+j*4+k*2+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  return m_data[i*2+j];
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
inline const X& tensor2d<X>::operator()(size_t i, size_t j) const
{
  if (i == j) return m_data[i];
  else        return m_zero[0];
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

template<class X> inline const X* tensor4 <X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor4 <X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor4 <X>::end()   const { return &m_data[0] + 16; }
template<class X> inline const X* tensor2 <X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2 <X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2 <X>::end()   const { return &m_data[0] + 4;  }
template<class X> inline const X* tensor2s<X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2s<X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2s<X>::end()   const { return &m_data[0] + 3;  }
template<class X> inline const X* tensor2d<X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2d<X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2d<X>::end()   const { return &m_data[0] + 2;  }
template<class X> inline const X* vector  <X>::data()  const { return &m_data[0];      }
template<class X> inline auto     vector  <X>::begin() const { return &m_data[0];      }
template<class X> inline auto     vector  <X>::end()   const { return &m_data[0] + 2;  }

// =================================================================================================
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
// tensor products
// =================================================================================================

template<class X>
inline reg::tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

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
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

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
inline reg::tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

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
inline reg::tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

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
inline reg::tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  reg::tensor2d<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  reg::vector<X> C(static_cast<X>(0));

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
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      for ( size_t l = 0 ; l < 2 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      for ( size_t l = 0 ; l < 2 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  reg::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t k = 0 ; k < 2 ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  reg::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

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

// -------------------------------------------------------------------------------------------------

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

// -------------------------------------------------------------------------------------------------

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

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor4<X>::RT() const
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor4<X> tensor4<X>::LT() const
{
  reg::tensor4<X> C;

  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      for ( size_t k = 0 ; k < 2 ; ++k )
        for ( size_t l = 0 ; l < 2 ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2<X> tensor2<X>::T() const
{
  reg::tensor2<X> C;

  C[0] = m_data[0];
  C[2] = m_data[1];
  C[1] = m_data[2];
  C[3] = m_data[3];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2s<X> tensor2s<X>::T() const
{
  reg::tensor2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::T() const
{
  reg::tensor2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];

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
inline reg::tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  reg::tensor2<X> C;

  C[0] =                      m_data[3] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] = static_cast<X>(-1) * m_data[2] / D;
  C[3] =                      m_data[0] / D;

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

  C[0] =                      m_data[2] / D;
  C[1] = static_cast<X>(-1) * m_data[1] / D;
  C[2] =                      m_data[0] / D;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline reg::tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  reg::tensor2d<X> C;

  C[0] = static_cast<X>(1) / m_data[0];
  C[1] = static_cast<X>(1) / m_data[1];

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
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor4<X>::operator== (const reg::tensor4<X> &B) const
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
inline bool tensor2<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;
  if ( m_data[3] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const reg::tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[1] ) return false;
  if ( m_data[3] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const reg::tensor2d<X> &B) const
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
inline bool tensor2s<X>::operator== (const reg::tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[1] != B[2] ) return false;
  if ( m_data[2] != B[3] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const reg::tensor2d<X> &B) const
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
inline bool tensor2d<X>::operator== (const reg::tensor2d<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const reg::tensor2<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[3] ) return false;
  if (              B[1] ) return false;
  if (              B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const reg::tensor2s<X> &B) const
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

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool vector<X>::operator== (const reg::vector<X> &B) const
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

}}} // namespace ...

#endif

