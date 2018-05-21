/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR4_CPP
#define CPPMAT_VAR_CARTESIAN_TENSOR4_CPP

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
tensor4<X>::tensor4(size_t nd) : cppmat::array<X>({nd,nd,nd,nd})
{
  ND = nd;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X>::tensor4(const cppmat::array<X> &A) : cppmat::array<X>(A)
{
  assert( this->mRank == 4 );

  ND = this->mShape[0];
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
tensor4<X> tensor4<X>::Random(size_t nd, X lower, X upper)
{
  tensor4<X> out(nd);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Arange(size_t nd)
{
  tensor4<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Zero(size_t nd)
{
  tensor4<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Ones(size_t nd)
{
  tensor4<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Constant(size_t nd, X D)
{
  tensor4<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Copy(size_t nd, const std::vector<X> &D)
{
  tensor4<X> out(nd);

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor4<X> tensor4<X>::Copy(size_t nd, Iterator first)
{
  tensor4<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
tensor4<X> tensor4<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  tensor4<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

template<class X>
inline
tensor4<X> tensor4<X>::I(size_t nd)
{
  tensor4<X> out(nd);

  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Irt(size_t nd)
{
  tensor4<X> out(nd);

  out.setIrt();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Id(size_t nd)
{
  tensor4<X> out(nd);

  out.setId();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Is(size_t nd)
{
  tensor4<X> out(nd);

  out.setIs();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::Isd(size_t nd)
{
  tensor4<X> out(nd);

  out.setIsd();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::II(size_t nd)
{
  tensor4<X> out(nd);

  out.setII();

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void tensor4<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::array<X>::resize({nd,nd,nd,nd});
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X>
inline
size_t tensor4<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize: identity tensors
// =================================================================================================

template<class X>
inline
void tensor4<X>::setI()
{
  cppmat::array<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          if ( i == l and j == k )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void tensor4<X>::setIrt()
{
  cppmat::array<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          if ( i == k and j == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void tensor4<X>::setIs()
{
  (*this) = ( tensor4<X>::I(ND) + tensor4<X>::Irt(ND) ) / static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void tensor4<X>::setId()
{
  (*this) = tensor4<X>::I(ND) - tensor4<X>::II(ND)/static_cast<X>(ND);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void tensor4<X>::setIsd()
{
  (*this) = tensor4<X>::Is(ND) - tensor4<X>::II(ND)/static_cast<X>(ND);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void tensor4<X>::setII()
{
  cppmat::array<X>::setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          if ( i == j and k == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
tensor4<X> tensor4<X>::ddot(const tensor4<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor4<X>::ddot(const tensor2<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor4<X>::ddot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor4<X>::ddot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
tensor4<X> tensor4<X>::T() const
{
  return cppmat::cartesian::T(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::RT() const
{
  return cppmat::cartesian::RT(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor4<X>::LT() const
{
  return cppmat::cartesian::LT(*this);
}

// =================================================================================================

}} // namespace ...

#endif

