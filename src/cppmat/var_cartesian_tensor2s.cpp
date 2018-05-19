/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */


#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2S_CPP
#define CPPMAT_VAR_CARTESIAN_TENSOR2S_CPP

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
tensor2s<X>::tensor2s(size_t nd) : cppmat::symmetric::matrix<X>(nd,nd)
{
  ND = nd;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X>::tensor2s(const cppmat::symmetric::matrix<X> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// -------------------------------------------------------------------------------------------------

#ifndef CPPMAT_NOCONVERT
template<class X>
inline
tensor2s<X>::tensor2s(const cppmat::diagonal::matrix<X> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}
#endif

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::Arange(size_t nd)
{
  tensor2s<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::Zero(size_t nd)
{
  tensor2s<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::Ones(size_t nd)
{
  tensor2s<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::Constant(size_t nd, X D)
{
  tensor2s<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::Copy(size_t nd, Iterator first)
{
  tensor2s<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  tensor2s<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::CopyDense(size_t nd, Iterator first)
{
  tensor2s<X> out(nd);

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::CopyDense(size_t nd, Iterator first, Iterator last)
{
  tensor2s<X> out(nd);

  out.setCopyDense(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::I(size_t nd)
{
  tensor2s<X> out(nd);

  out.setI();

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void tensor2s<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::symmetric::matrix<X>::resize(nd,nd);
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X>
inline
size_t tensor2s<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void tensor2s<X>::setI()
{
  cppmat::symmetric::matrix<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
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
X tensor2s<X>::ddot(const tensor2<X> &B) const
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
X tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += ( this->mData[i*ND-(i-1)*i/2] * B[i*ND-(i-1)*i/2] );

  for ( size_t i = 0 ; i<ND ; ++i )
    for ( size_t j = i+1 ; j<ND ; ++j )
      C += ( static_cast<X>(2) * this->mData[i*ND-(i-1)*i/2+j-i] * B[i*ND-(i-1)*i/2+j-i] );

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += this->mData[i*ND-(i-1)*i/2]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
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
tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
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
tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
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
vector<X> tensor2s<X>::dot(const vector<X> &B) const
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
tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
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
tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
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
tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
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
tensor2s<X> tensor2s<X>::T() const
{
  tensor2s<X> C = (*this);

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X tensor2s<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += (*this)(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::det() const
{
  if ( ND==2 )
   return this->mData[0] * this->mData[2] - this->mData[1] * this->mData[1];

  if ( ND==3 )
    return (                     this->mData[0] * this->mData[3] * this->mData[5] +
             static_cast<X>(2) * this->mData[1] * this->mData[2] * this->mData[4] ) -
           (                     this->mData[4] * this->mData[4] * this->mData[0] +
                                 this->mData[2] * this->mData[2] * this->mData[3] +
                                 this->mData[1] * this->mData[1] * this->mData[5] );

  throw std::runtime_error("'det' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::inv() const
{
  // compute determinant
  X D = det();

  // allocate result
  tensor2s<X> C(ND);

  if ( ND==2 )
  {
    C[0] =                      this->mData[2] / D;
    C[1] = static_cast<X>(-1) * this->mData[1] / D;
    C[2] =                      this->mData[0] / D;

    return C;
  }

  if ( ND==3 )
  {
    C[0] = (this->mData[3]*this->mData[5]-this->mData[4]*this->mData[4]) / D;
    C[1] = (this->mData[2]*this->mData[4]-this->mData[1]*this->mData[5]) / D;
    C[2] = (this->mData[1]*this->mData[4]-this->mData[2]*this->mData[3]) / D;
    C[3] = (this->mData[0]*this->mData[5]-this->mData[2]*this->mData[2]) / D;
    C[4] = (this->mData[2]*this->mData[1]-this->mData[0]*this->mData[4]) / D;
    C[5] = (this->mData[0]*this->mData[3]-this->mData[1]*this->mData[1]) / D;

    return C;
  }

  throw std::runtime_error("'inv' only implemented in 2D/3D, use e.g. 'Eigen'");
}

// =================================================================================================

}} // namespace ...

#endif

