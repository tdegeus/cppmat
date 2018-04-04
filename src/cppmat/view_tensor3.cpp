/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VIEW_TENSOR3_CPP
#define CPPMAT_VIEW_TENSOR3_CPP

#include "view_tensor3.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian3d {

// -------------------------------------------------------------------------------------------------

// alias
namespace normal = cppmat::cartesian3d;

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline tensor4<X>::tensor4()
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor4<X>::tensor4(const X *D)
{
  m_data = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>::tensor2()
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>::tensor2(const X *D)
{
  m_data = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>::tensor2s()
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>::tensor2s(const X *D)
{
  m_data = D;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d()
{
  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>::tensor2d(const X *D)
{
  m_data = D;

  m_zero[0] = static_cast<X>(0);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector()
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>::vector(const X *D)
{
  m_data = D;
}

// =================================================================================================
// assignment operators
// =================================================================================================

template<class X>
inline tensor4<X>& tensor4<X>::operator= (const tensor4<X> &D)
{
  m_data = &D[0];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2<X>& tensor2<X>::operator= (const tensor2<X> &D)
{
  m_data = &D[0];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2s<X>& tensor2s<X>::operator= (const tensor2s<X> &D)
{
  m_data = &D[0];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline tensor2d<X>& tensor2d<X>::operator= (const tensor2d<X> &D)
{
  m_data = &D[0];

  return *this;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline vector<X>& vector<X>::operator= (const vector<X> &D)
{
  m_data = &D[0];

  return *this;
}

// =================================================================================================
// map external pointer
// =================================================================================================

template<class X> inline void tensor4 <X>::map(const X *D) { m_data = D; }
template<class X> inline void tensor2 <X>::map(const X *D) { m_data = D; }
template<class X> inline void tensor2s<X>::map(const X *D) { m_data = D; }
template<class X> inline void tensor2d<X>::map(const X *D) { m_data = D; }
template<class X> inline void vector  <X>::map(const X *D) { m_data = D; }

// =================================================================================================
// cast to different class / type
// =================================================================================================

template<>
template<>
inline normal::tensor2s<double> tensor2<double>::cast<normal::tensor2s<double>>() const
{
  normal::tensor2s<double> out;

  out[0] =   m_data[0];
  out[1] = ( m_data[1] + m_data[3] ) / 2.;
  out[2] = ( m_data[2] + m_data[6] ) / 2.;
  out[3] =   m_data[4];
  out[4] = ( m_data[5] + m_data[7] ) / 2.;
  out[5] =   m_data[8];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline normal::tensor2d<double> tensor2<double>::cast<normal::tensor2d<double>>() const
{
  normal::tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[4];
  out[2] = m_data[8];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline normal::tensor2<double> tensor2s<double>::cast<normal::tensor2<double>>() const
{
  normal::tensor2<double> out;

  out[0]          = m_data[0];
  out[1] = out[3] = m_data[1];
  out[2] = out[6] = m_data[2];
  out[4]          = m_data[3];
  out[5] = out[7] = m_data[4];
  out[8]          = m_data[5];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline normal::tensor2d<double> tensor2s<double>::cast<normal::tensor2d<double>>() const
{
  normal::tensor2d<double> out;

  out[0] = m_data[0];
  out[1] = m_data[3];
  out[2] = m_data[5];

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline normal::tensor2<double> tensor2d<double>::cast<normal::tensor2<double>>() const
{
  normal::tensor2<double> out;

  out[0] = m_data[0];
  out[4] = m_data[1];
  out[8] = m_data[2];

  out[1] = out[2] = out[3] = out[5] = out[6] = out[7] = 0.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

template<>
template<>
inline normal::tensor2s<double> tensor2d<double>::cast<normal::tensor2s<double>>() const
{
  normal::tensor2s<double> out;

  out[0] = m_data[0];
  out[3] = m_data[1];
  out[5] = m_data[2];

  out[1] = out[2] = out[4] = 0.0;

  return out;
}

// =================================================================================================
// get dimensions
// =================================================================================================

template<class X> inline size_t tensor4 <X>::size() const { return 81; }
template<class X> inline size_t tensor2 <X>::size() const { return 9;  }
template<class X> inline size_t tensor2s<X>::size() const { return 6;  }
template<class X> inline size_t tensor2d<X>::size() const { return 3;  }
template<class X> inline size_t vector  <X>::size() const { return 3;  }
template<class X> inline size_t tensor4 <X>::ndim() const { return 3;  }
template<class X> inline size_t tensor2 <X>::ndim() const { return 3;  }
template<class X> inline size_t tensor2s<X>::ndim() const { return 3;  }
template<class X> inline size_t tensor2d<X>::ndim() const { return 3;  }
template<class X> inline size_t vector  <X>::ndim() const { return 3;  }

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::shape() const
{
  std::vector<size_t> out(4,3);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor2<X>::shape() const
{
  std::vector<size_t> out(2,3);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> vector<X>::shape() const
{
  std::vector<size_t> out(1,3);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::vector<size_t> tensor4<X>::strides(bool bytes) const
{
  std::vector<size_t> out = { 27, 9, 3, 1 };

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
  std::vector<size_t> out = { 3, 1 };

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
  return m_data[i*27+j*9+k*3+l];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2<X>::operator()(size_t i, size_t j) const
{
  return m_data[i*3+j];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline const X& tensor2s<X>::operator()(size_t i, size_t j) const
{
  if ( i == 0 ) {
    if ( j == 0 ) return m_data[0];
    if ( j == 1 ) return m_data[1];
    else          return m_data[2];
  }
  if ( i == 1 ) {
    if ( j == 0 ) return m_data[1];
    if ( j == 1 ) return m_data[3];
    else          return m_data[4];
  }
  else {
    if ( j == 0 ) return m_data[2];
    if ( j == 1 ) return m_data[4];
    else          return m_data[5];
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
template<class X> inline auto     tensor4 <X>::end()   const { return &m_data[0] + 81; }
template<class X> inline const X* tensor2 <X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2 <X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2 <X>::end()   const { return &m_data[0] + 9;  }
template<class X> inline const X* tensor2s<X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2s<X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2s<X>::end()   const { return &m_data[0] + 6;  }
template<class X> inline const X* tensor2d<X>::data()  const { return &m_data[0];      }
template<class X> inline auto     tensor2d<X>::begin() const { return &m_data[0];      }
template<class X> inline auto     tensor2d<X>::end()   const { return &m_data[0] + 3;  }
template<class X> inline const X* vector  <X>::data()  const { return &m_data[0];      }
template<class X> inline auto     vector  <X>::begin() const { return &m_data[0];      }
template<class X> inline auto     vector  <X>::end()   const { return &m_data[0] + 3;  }

// =================================================================================================
// arithmetic operators
// =================================================================================================

template<class X>
inline normal::tensor4<X> operator* (const tensor4<X> &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator/ (const tensor4<X> &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator+ (const tensor4<X> &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator- (const tensor4<X> &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator* (const tensor4<X> &A, const X &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator/ (const tensor4<X> &A, const X &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator+ (const tensor4<X> &A, const X &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator- (const tensor4<X> &A, const X &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator* (const X &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator/ (const X &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator+ (const X &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> operator- (const X &A, const tensor4<X> &B)
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator* (const tensor2<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] * B[i  ];
    C[i+1] = A[i+1] * B[i+1];
    C[i+2] = A[i+2] * B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator/ (const tensor2<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] / B[i  ];
    C[i+1] = A[i+1] / B[i+1];
    C[i+2] = A[i+2] / B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] + B[i  ];
    C[i+1] = A[i+1] + B[i+1];
    C[i+2] = A[i+2] + B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] - B[i  ];
    C[i+1] = A[i+1] - B[i+1];
    C[i+2] = A[i+2] - B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator* (const tensor2<X> &A, const tensor2s<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[3] * B[1];
  C[2] = A[2] * B[2]; C[6] = A[6] * B[2];
  C[4] = A[4] * B[3];
  C[5] = A[5] * B[4]; C[7] = A[7] * B[4];
  C[8] = A[8] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator/ (const tensor2<X> &A, const tensor2s<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[3] / B[1];
  C[2] = A[2] / B[2]; C[6] = A[6] / B[2];
  C[4] = A[4] / B[3];
  C[5] = A[5] / B[4]; C[7] = A[7] / B[4];
  C[8] = A[8] / B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2<X> &A, const tensor2s<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[3] + B[1];
  C[2] = A[2] + B[2]; C[6] = A[6] + B[2];
  C[4] = A[4] + B[3];
  C[5] = A[5] + B[4]; C[7] = A[7] + B[4];
  C[8] = A[8] + B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2<X> &A, const tensor2s<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[3] - B[1];
  C[2] = A[2] - B[2]; C[6] = A[6] - B[2];
  C[4] = A[4] - B[3];
  C[5] = A[5] - B[4]; C[7] = A[7] - B[4];
  C[8] = A[8] - B[5];

  return C;
}


// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2<X> &A, const tensor2d<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3];
  C[4] = A[4] + B[1];
  C[5] = A[5];
  C[6] = A[6];
  C[7] = A[7];
  C[8] = A[8] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2<X> &A, const tensor2d<X> &B)
{
  normal::tensor2 <X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3];
  C[4] = A[4] - B[1];
  C[5] = A[5];
  C[6] = A[6];
  C[7] = A[7];
  C[8] = A[8] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator* (const tensor2<X> &A, const X &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] * B;
    C[i+1] = A[i+1] * B;
    C[i+2] = A[i+2] * B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator/ (const tensor2<X> &A, const X &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] / B;
    C[i+1] = A[i+1] / B;
    C[i+2] = A[i+2] / B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2<X> &A, const X &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] + B;
    C[i+1] = A[i+1] + B;
    C[i+2] = A[i+2] + B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2<X> &A, const X &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A[i  ] - B;
    C[i+1] = A[i+1] - B;
    C[i+2] = A[i+2] - B;
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator* (const tensor2s<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1]; C[3] = A[1] * B[3];
  C[2] = A[2] * B[2]; C[6] = A[2] * B[6];
  C[4] = A[3] * B[4];
  C[5] = A[4] * B[5]; C[7] = A[4] * B[7];
  C[8] = A[5] * B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator/ (const tensor2s<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1]; C[3] = A[1] / B[3];
  C[2] = A[2] / B[2]; C[6] = A[2] / B[6];
  C[4] = A[3] / B[4];
  C[5] = A[4] / B[5]; C[7] = A[4] / B[7];
  C[8] = A[5] / B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2s<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1]; C[3] = A[1] + B[3];
  C[2] = A[2] + B[2]; C[6] = A[2] + B[6];
  C[4] = A[3] + B[4];
  C[5] = A[4] + B[5]; C[7] = A[4] + B[7];
  C[8] = A[5] + B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2s<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1]; C[3] = A[1] - B[3];
  C[2] = A[2] - B[2]; C[6] = A[2] - B[6];
  C[4] = A[3] - B[4];
  C[5] = A[4] - B[5]; C[7] = A[4] - B[7];
  C[8] = A[5] - B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const tensor2d<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] =        B[3];
  C[4] = A[1] + B[4];
  C[5] =        B[5];
  C[6] =        B[6];
  C[7] =        B[7];
  C[8] = A[2] + B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const tensor2d<X> &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] =      - B[3];
  C[4] = A[1] - B[4];
  C[5] =      - B[5];
  C[6] =      - B[6];
  C[7] =      - B[7];
  C[8] = A[2] - B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator* (const X &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A * B[i  ];
    C[i+1] = A * B[i+1];
    C[i+2] = A * B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator/ (const X &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A / B[i  ];
    C[i+1] = A / B[i+1];
    C[i+2] = A / B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator+ (const X &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A + B[i  ];
    C[i+1] = A + B[i+1];
    C[i+2] = A + B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> operator- (const X &A, const tensor2<X> &B)
{
  normal::tensor2<X> C;

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C[i  ] = A - B[i  ];
    C[i+1] = A - B[i+1];
    C[i+2] = A - B[i+2];
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator* (const tensor2s<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];
  C[3] = A[3] * B[3];
  C[4] = A[4] * B[4];
  C[5] = A[5] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator/ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];
  C[3] = A[3] / B[3];
  C[4] = A[4] / B[4];
  C[5] = A[5] / B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];
  C[3] = A[3] + B[3];
  C[4] = A[4] + B[4];
  C[5] = A[5] + B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];
  C[3] = A[3] - B[3];
  C[4] = A[4] - B[4];
  C[5] = A[5] - B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const tensor2s<X> &A, const tensor2d<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] + B[1];
  C[4] = A[4];
  C[5] = A[5] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const tensor2s<X> &A, const tensor2d<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1];
  C[2] = A[2];
  C[3] = A[3] - B[1];
  C[4] = A[4];
  C[5] = A[5] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator* (const tensor2s<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;
  C[3] = A[3] * B;
  C[4] = A[4] * B;
  C[5] = A[5] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator/ (const tensor2s<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;
  C[3] = A[3] / B;
  C[4] = A[4] / B;
  C[5] = A[5] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const tensor2s<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;
  C[3] = A[3] + B;
  C[4] = A[4] + B;
  C[5] = A[5] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const tensor2s<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;
  C[3] = A[3] - B;
  C[4] = A[4] - B;
  C[5] = A[5] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const tensor2d<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] + B;
  C[1] =        B;
  C[2] =        B;
  C[3] = A[1] + B;
  C[4] =        B;
  C[5] = A[2] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const tensor2d<X> &A, const X &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] - B;
  C[1] =      - B;
  C[2] =      - B;
  C[3] = A[1] - B;
  C[4] =      - B;
  C[5] = A[2] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] + B[0];
  C[1] =        B[1];
  C[2] =        B[2];
  C[3] = A[1] + B[3];
  C[4] =        B[4];
  C[5] = A[2] + B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const tensor2d<X> &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A[0] - B[0];
  C[1] =      - B[1];
  C[2] =      - B[2];
  C[3] = A[1] - B[3];
  C[4] =      - B[4];
  C[5] = A[2] - B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator* (const X &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];
  C[3] = A * B[3];
  C[4] = A * B[4];
  C[5] = A * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator/ (const X &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];
  C[3] = A / B[3];
  C[4] = A / B[4];
  C[5] = A / B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const X &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];
  C[3] = A + B[3];
  C[4] = A + B[4];
  C[5] = A + B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const X &A, const tensor2s<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];
  C[3] = A - B[3];
  C[4] = A - B[4];
  C[5] = A - B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator+ (const X &A, const tensor2d<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A + B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A + B[1];
  C[4] = A;
  C[5] = A + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> operator- (const X &A, const tensor2d<X> &B)
{
  normal::tensor2s<X> C;

  C[0] = A - B[0];
  C[1] = A;
  C[2] = A;
  C[3] = A - B[1];
  C[4] = A;
  C[5] = A - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator+ (const tensor2d<X> &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator- (const tensor2d<X> &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2 <X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[4];
  C[2] = A[2] * B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2 <X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[4];
  C[2] = A[2] / B[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2d<X> &A, const tensor2s<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[3];
  C[2] = A[2] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator/ (const tensor2d<X> &A, const tensor2s<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[3];
  C[2] = A[2] / B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2d<X> &A, const X &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator/ (const tensor2d<X> &A, const X &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2 <X> &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[4] * B[1];
  C[2] = A[8] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const tensor2s<X> &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[3] * B[1];
  C[2] = A[5] * B[2];

  return C;
}


// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> operator* (const X &A, const tensor2d<X> &B)
{
  normal::tensor2d<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator* (const vector<X> &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A[0] * B[0];
  C[1] = A[1] * B[1];
  C[2] = A[2] * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator/ (const vector<X> &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A[0] / B[0];
  C[1] = A[1] / B[1];
  C[2] = A[2] / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator+ (const vector<X> &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator- (const vector<X> &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A[0] - B[0];
  C[1] = A[1] - B[1];
  C[2] = A[2] - B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator* (const vector<X> &A, const X &B)
{
  normal::vector<X> C;

  C[0] = A[0] * B;
  C[1] = A[1] * B;
  C[2] = A[2] * B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator/ (const vector<X> &A, const X &B)
{
  normal::vector<X> C;

  C[0] = A[0] / B;
  C[1] = A[1] / B;
  C[2] = A[2] / B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator+ (const vector<X> &A, const X &B)
{
  normal::vector<X> C;

  C[0] = A[0] + B;
  C[1] = A[1] + B;
  C[2] = A[2] + B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator- (const vector<X> &A, const X &B)
{
  normal::vector<X> C;

  C[0] = A[0] - B;
  C[1] = A[1] - B;
  C[2] = A[2] - B;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator* (const X &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A * B[0];
  C[1] = A * B[1];
  C[2] = A * B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator/ (const X &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A / B[0];
  C[1] = A / B[1];
  C[2] = A / B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator+ (const X &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A + B[0];
  C[1] = A + B[1];
  C[2] = A + B[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> operator- (const X &A, const vector<X> &B)
{
  normal::vector<X> C;

  C[0] = A - B[0];
  C[1] = A - B[1];
  C[2] = A - B[2];

  return C;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline normal::tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          for ( size_t m = 0 ; m < 3 ; ++m )
            for ( size_t n = 0 ; n < 3 ; ++n )
              C(i,j,m,n) += (*this)(i,j,k,l) * B(l,k,m,n);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j) += (*this)(i,j,k,l) * B(l,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,j) += (*this)(i,j,k,k) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[4]*B[1] + m_data[8]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
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
  C += m_data[2] * B[2] * static_cast<X>(2);
  C += m_data[3] * B[3];
  C += m_data[4] * B[4] * static_cast<X>(2);
  C += m_data[5] * B[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[3]*B[1] + m_data[5]*B[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      for ( size_t l = 0 ; l < 3 ; ++l )
        C(k,l) += m_data[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[4] + m_data[2]*B[8];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[3] + m_data[2]*B[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  return m_data[0]*B[0] + m_data[1]*B[1] + m_data[2]*B[2];
}


// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  normal::tensor2d<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  normal::vector<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    C(i) += (*this)(i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::dot(const vector<X> &B) const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 3 ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      for ( size_t l = 0 ; l < 3 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      for ( size_t l = 0 ; l < 3 ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  normal::tensor4<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t k = 0 ; k < 3 ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  normal::tensor2<X> C(static_cast<X>(0));

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> vector<X>::cross(const vector<X> &B) const
{
  normal::vector<X> C;

  C[0] =                     m_data[1]*B[2]-B[1]*m_data[2] ;
  C[1] = static_cast<X>(-1)*(m_data[0]*B[2]-B[0]*m_data[2]);
  C[2] =                     m_data[0]*B[1]-B[0]*m_data[1] ;

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> ddot(const tensor4<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor4<X> &A, const tensor2<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor4<X> &A, const tensor2s<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor4<X> &A, const tensor2d<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor2<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor2s<X> &A, const tensor4<X> &B)
{
  return A.ddot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> ddot(const tensor2d<X> &A, const tensor4<X> &B)
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
inline normal::tensor2<X> dot(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dot(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> dot(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const tensor2<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const tensor2s<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const tensor2d<X> &A, const vector<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const vector<X> &A, const tensor2<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const vector<X> &A, const tensor2s<X> &B)
{
  return A.dot(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> dot(const vector<X> &A, const tensor2d<X> &B)
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
inline normal::tensor4<X> dyadic(const tensor2<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2s<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2s<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> dyadic(const tensor2d<X> &A, const tensor2d<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> dyadic(const vector<X> &A, const vector<X> &B)
{
  return A.dyadic(B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::vector<X> cross(const vector<X> &A, const vector<X> &B)
{
  return A.cross (B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline normal::tensor4<X> tensor4<X>::T() const
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(l,k,j,i) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor4<X>::RT() const
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(i,j,l,k) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> tensor4<X>::LT() const
{
  normal::tensor4<X> C;

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          C(j,i,k,l) = (*this)(i,j,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::T() const
{
  normal::tensor2<X> C;

  C[0] = m_data[0];
  C[3] = m_data[1];
  C[6] = m_data[2];
  C[1] = m_data[3];
  C[4] = m_data[4];
  C[7] = m_data[5];
  C[2] = m_data[6];
  C[5] = m_data[7];
  C[8] = m_data[8];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> tensor2s<X>::T() const
{
  normal::tensor2s<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];
  C[3] = m_data[3];
  C[4] = m_data[4];
  C[5] = m_data[5];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> tensor2d<X>::T() const
{
  normal::tensor2d<X> C;

  C[0] = m_data[0];
  C[1] = m_data[1];
  C[2] = m_data[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> transpose(const tensor2<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> transpose(const tensor2s<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> transpose(const tensor2d<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> transpose(const tensor4<X> &A)
{
  return A.T();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> transposeR(const tensor4<X> &A)
{
  return A.RT();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor4<X> transposeL(const tensor4<X> &A)
{
  return A.LT();
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline X tensor2<X>::trace() const
{
  return m_data[0] + m_data[4] + m_data[8];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::trace() const
{
  return m_data[0] + m_data[3] + m_data[5];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::trace() const
{
  return m_data[0] + m_data[1] + m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::det() const
{
  return ( m_data[0] * m_data[4] * m_data[8] +
           m_data[1] * m_data[5] * m_data[6] +
           m_data[2] * m_data[3] * m_data[7] ) -
         ( m_data[2] * m_data[4] * m_data[6] +
           m_data[1] * m_data[3] * m_data[8] +
           m_data[0] * m_data[5] * m_data[7] );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2s<X>::det() const
{
  return (                     m_data[0] * m_data[3] * m_data[5] +
           static_cast<X>(2) * m_data[1] * m_data[2] * m_data[4] ) -
         (                     m_data[4] * m_data[4] * m_data[0] +
                               m_data[2] * m_data[2] * m_data[3] +
                               m_data[1] * m_data[1] * m_data[5] );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::det() const
{
  return m_data[0] * m_data[1] * m_data[2];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  normal::tensor2<X> C;

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

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  normal::tensor2s<X> C;

  C[0] = (m_data[3]*m_data[5]-m_data[4]*m_data[4]) / D;
  C[1] = (m_data[2]*m_data[4]-m_data[1]*m_data[5]) / D;
  C[2] = (m_data[1]*m_data[4]-m_data[2]*m_data[3]) / D;
  C[3] = (m_data[0]*m_data[5]-m_data[2]*m_data[2]) / D;
  C[4] = (m_data[2]*m_data[1]-m_data[0]*m_data[4]) / D;
  C[5] = (m_data[0]*m_data[3]-m_data[1]*m_data[1]) / D;
  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  normal::tensor2d<X> C;

  C[0] = static_cast<X>(1) / m_data[0];
  C[1] = static_cast<X>(1) / m_data[1];
  C[2] = static_cast<X>(1) / m_data[2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2<X> inv(const tensor2<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2s<X> inv(const tensor2s<X> &A)
{
  return A.inv();
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline normal::tensor2d<X> inv(const tensor2d<X> &A)
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
  for ( size_t i = 0 ; i < 81 ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2<X> &B) const
{
  for ( size_t i = 0 ; i < 9 ; ++i )
    if ( m_data[i] != B[i] )
      return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2s<X> &B) const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( m_data[i*3+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::operator== (const tensor2d<X> &B) const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( m_data[i*3+j] != B(i,j) )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2s<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;
  if ( m_data[3] != B[3] ) return false;
  if ( m_data[4] != B[4] ) return false;
  if ( m_data[5] != B[5] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2<X> &B) const
{
  for ( size_t i = 0 ; i < 3 ; ++i ) {
    for ( size_t j = i ; j < 3 ; ++j ) {
      if ( m_data[i*3-(i-1)*i/2+j-i] != B(i,j) ) return false;
      if ( m_data[i*3-(i-1)*i/2+j-i] != B(j,i) ) return false;
    }
  }

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::operator== (const tensor2d<X> &B) const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i ; j < 3 ; ++j )
      if ( m_data[i*3-(i-1)*i/2+j-i] != B(i,j) ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2d<X> &B) const
{
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2d<X>::operator== (const tensor2<X> &B) const
{
  for ( size_t i = 0 ; i < 3 ; ++i ) {
    for ( size_t j = 0 ; j < 3 ; ++j ) {
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
  for ( size_t i = 0 ; i < 3 ; ++i ) {
    for ( size_t j = i ; j < 3 ; ++j ) {
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
  if ( m_data[0] != B[0] ) return false;
  if ( m_data[1] != B[1] ) return false;
  if ( m_data[2] != B[2] ) return false;
  return true;
}

// =================================================================================================
// structure check
// =================================================================================================

template<class X>
inline bool tensor2<X>::issymmetric() const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      if ( m_data[i*3+j] != m_data[j*3+i] )
        return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      if ( i != j )
        if ( m_data[i*3+j] )
          return false;

  return true;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline bool tensor2s<X>::isdiagonal() const
{
  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = i+1 ; j < 3 ; ++j )
      if ( m_data[i*3-(i-1)*i/2+j-i] )
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

  for ( size_t i = 0 ; i < 81 ; i += 3 )
  {
    C += std::abs(m_data[i  ]);
    C += std::abs(m_data[i+1]);
    C += std::abs(m_data[i+2]);
  }

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2<X>::norm() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < 9 ; i += 3 )
  {
    C += std::abs(m_data[i  ]);
    C += std::abs(m_data[i+1]);
    C += std::abs(m_data[i+2]);
  }

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
  C += std::abs(m_data[3]);
  C += std::abs(m_data[4]);
  C += std::abs(m_data[5]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X tensor2d<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);
  C += std::abs(m_data[2]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::norm() const
{
  X C;

  C  = std::abs(m_data[0]);
  C += std::abs(m_data[1]);
  C += std::abs(m_data[2]);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline X vector<X>::length() const
{
  X C;

  C  = std::pow(m_data[0],2.);
  C += std::pow(m_data[1],2.);
  C += std::pow(m_data[2],2.);

  return std::sqrt(C);
}

// =================================================================================================
// formatted print
// =================================================================================================

template<class X>
inline void tensor4<X>::printf(std::string fmt) const
{
  std::string gmt = std::to_string(std::to_string(3).size());
  fmt = "(%"+gmt+"d,%"+gmt+"d,%"+gmt+"d,%"+gmt+"d) "+fmt+"\n";

  for ( size_t i = 0 ; i < 3 ; ++i )
    for ( size_t j = 0 ; j < 3 ; ++j )
      for ( size_t k = 0 ; k < 3 ; ++k )
        for ( size_t l = 0 ; l < 3 ; ++l )
          std::printf(fmt.c_str(), i, j, k, l, (*this)(i,j,k,l) );
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(0,0),(*this)(0,1),(*this)(0,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(1,0),(*this)(1,1),(*this)(1,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(2,0),(*this)(2,1),(*this)(2,2));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2s<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(0,0),(*this)(0,1),(*this)(0,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(1,0),(*this)(1,1),(*this)(1,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(2,0),(*this)(2,1),(*this)(2,2));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void tensor2d<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(0,0),(*this)(0,1),(*this)(0,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(1,0),(*this)(1,1),(*this)(1,2));
  std::printf((fmt + "," + fmt + "," + fmt + ";\n").c_str(),(*this)(2,0),(*this)(2,1),(*this)(2,2));
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline void vector<X>::printf(std::string fmt) const
{
  std::printf((fmt + "," + fmt + "," + fmt + "\n").c_str(), m_data[0], m_data[1], m_data[2]);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, tensor4<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  for ( size_t i = 0 ; i < 3 ; ++i ) {
    for ( size_t j = 0 ; j < 3 ; ++j ) {
      for ( size_t k = 0 ; k < 3 ; ++k ) {
        for ( size_t l = 0 ; l < 3 ; ++l ) {
          out << "(" << i << "," << j << "," << k << "," << l << ") = ";
          out << std::setw(w) << std::setprecision(p) << src(i,j,k,l);
          if ( !(i==2 and j==2 and k==2 and l==2) ) out << std::endl;
        }
      }
    }
  }

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, tensor2<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(2,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,2) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, tensor2s<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(2,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,2) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, tensor2d<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(0,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(1,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1,2) << ";" << std::endl;

  out << std::setw(w) << std::setprecision(p) << src(2,0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2,2) << ";";

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline std::ostream& operator<<(std::ostream& out, vector<X>& src)
{
  auto w = out.width();
  auto p = out.precision();

  out << std::setw(w) << std::setprecision(p) << src(0) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(1) << ", ";
  out << std::setw(w) << std::setprecision(p) << src(2);

  return out;
}

// =================================================================================================

}}} // namespace ...

#endif

