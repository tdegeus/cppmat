/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2_CPP
#define CPPMAT_VAR_CARTESIAN_TENSOR2_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// constructors
// =================================================================================================

template<class X>
inline
tensor2<X>::tensor2(size_t nd) : cppmat::matrix<X>(nd,nd)
{
  ND = nd;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X>::tensor2(const cppmat::array<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->M;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X>::tensor2(const cppmat::matrix<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->M;
}

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline
tensor2<X>::tensor2(const cppmat::symmetric::matrix<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->M;
}
#endif

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline
tensor2<X>::tensor2(const cppmat::diagonal::matrix<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->M;
}
#endif

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X>::tensor2(size_t nd, const std::vector<X> &D) : cppmat::matrix<X>(nd,nd, D)
{
  ND = nd;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::Arange(size_t nd)
{
  tensor2<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::Zero(size_t nd)
{
  tensor2<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::Ones(size_t nd)
{
  tensor2<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::Constant(size_t nd, X D)
{
  tensor2<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2<X> tensor2<X>::Copy(size_t nd, Iterator first)
{
  tensor2<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2<X> tensor2<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  tensor2<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::I(size_t nd)
{
  tensor2<X> out(nd);

  out.setI();

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void tensor2<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::matrix<X>::resize(nd,nd);
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X>
inline
size_t tensor2<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void tensor2<X>::setI()
{
  cppmat::matrix<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(k,l) += (*this)(i,j) * B(j,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C += (*this)(i,j) * B(j,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,k) += (*this)(i,j) * B(j,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += (*this)(i,j) * B(j,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i) += (*this)(i,j) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          C(i,j,k,l) += (*this)(i,j) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        C(i,j,k,k) += (*this)(i,j) * B(k,k);

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2<X>::T() const
{
  tensor2<X> C(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j,i) = (*this)(i,j);

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X tensor2<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::det() const
{
  if ( ND==2 )
   return this->mData[0] * this->mData[3] - this->mData[1] * this->mData[2];

  if ( ND==3 )
    return ( this->mData[0] * this->mData[4] * this->mData[8] +
             this->mData[1] * this->mData[5] * this->mData[6] +
             this->mData[2] * this->mData[3] * this->mData[7] ) -
           ( this->mData[2] * this->mData[4] * this->mData[6] +
             this->mData[1] * this->mData[3] * this->mData[8] +
             this->mData[0] * this->mData[5] * this->mData[7] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2<X> C(ND);

  if ( ND==2 )
  {
    C[0] =                      this->mData[3] / D;
    C[1] = static_cast<X>(-1) * this->mData[1] / D;
    C[2] = static_cast<X>(-1) * this->mData[2] / D;
    C[3] =                      this->mData[0] / D;

    return C;
  }

  if ( ND==3 )
  {
    C[0] = (this->mData[4]*this->mData[8]-this->mData[5]*this->mData[7]) / D;
    C[1] = (this->mData[2]*this->mData[7]-this->mData[1]*this->mData[8]) / D;
    C[2] = (this->mData[1]*this->mData[5]-this->mData[2]*this->mData[4]) / D;
    C[3] = (this->mData[5]*this->mData[6]-this->mData[3]*this->mData[8]) / D;
    C[4] = (this->mData[0]*this->mData[8]-this->mData[2]*this->mData[6]) / D;
    C[5] = (this->mData[2]*this->mData[3]-this->mData[0]*this->mData[5]) / D;
    C[6] = (this->mData[3]*this->mData[7]-this->mData[4]*this->mData[6]) / D;
    C[7] = (this->mData[1]*this->mData[6]-this->mData[0]*this->mData[7]) / D;
    C[8] = (this->mData[0]*this->mData[4]-this->mData[1]*this->mData[3]) / D;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// =================================================================================================

}} // namespace ...

#endif

