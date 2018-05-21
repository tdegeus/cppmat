/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2D_CPP
#define CPPMAT_VAR_CARTESIAN_TENSOR2D_CPP

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
tensor2d<X>::tensor2d() : cppmat::diagonal::matrix<X>()
{
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X>::tensor2d(size_t nd) : cppmat::diagonal::matrix<X>(nd,nd)
{
  ND = nd;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X>::tensor2d(const cppmat::diagonal::matrix<X> &A) : cppmat::diagonal::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
tensor2d<X> tensor2d<X>::Random(size_t nd, X lower, X upper)
{
  tensor2d<X> out(nd);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::Arange(size_t nd)
{
  tensor2d<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::Zero(size_t nd)
{
  tensor2d<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::Ones(size_t nd)
{
  tensor2d<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::Constant(size_t nd, X D)
{
  tensor2d<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2d<X> tensor2d<X>::Copy(size_t nd, Iterator first)
{
  tensor2d<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2d<X> tensor2d<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  tensor2d<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2d<X> tensor2d<X>::CopyDense(size_t nd, Iterator first)
{
  tensor2d<X> out(nd);

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor2d<X> tensor2d<X>::CopyDense(size_t nd, Iterator first, Iterator last)
{
  tensor2d<X> out(nd);

  out.setCopyDense(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::I(size_t nd)
{
  tensor2d<X> out(nd);

  out.setI();

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void tensor2d<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::diagonal::matrix<X>::resize(nd,nd);
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X>
inline
size_t tensor2d<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize
// =================================================================================================

template<class X>
inline
void tensor2d<X>::setI()
{
  cppmat::diagonal::matrix<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    this->mData[i] = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2d<X>::ddot(const tensor4<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(k,l) += this->mData[i]*B(i,i,k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2d<X>::ddot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += this->mData[i]*B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2d<X>::ddot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += this->mData[i]*B[i*ND-(i-1)*i/2];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2d<X>::ddot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += this->mData[i]*B[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2d<X>::dot(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2d<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,k) += (*this)(i,i) * B(i,k);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2d<X> C = tensor2d<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i,i) += (*this)(i,i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> tensor2d<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += (*this)(i,i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2d<X>::dyadic(const tensor2<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2d<X>::dyadic(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      for ( size_t l = 0 ; l < ND ; ++l )
        C(i,i,k,l) += (*this)(i,i) * B(k,l);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2d<X>::dyadic(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor4<X> C = tensor4<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t k = 0 ; k < ND ; ++k )
      C(i,i,k,k) += (*this)(i,i) * B(k,k);

  return C;
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
tensor2d<X> tensor2d<X>::T() const
{
  tensor2d<X> C = (*this);

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X tensor2d<X>::trace() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2d<X>::det() const
{
  X C = static_cast<X>(1);

  for ( size_t i = 0 ; i < ND ; ++i )
    C *= (*this)[i];

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2d<X> tensor2d<X>::inv() const
{
  // allocate result
  tensor2d<X> C(ND);

  for ( size_t i = 0; i < ND ; ++i )
    C[i] = static_cast<X>(1) / this->mData[i];

  return C;
}


// =================================================================================================

}} // namespace ...

#endif

