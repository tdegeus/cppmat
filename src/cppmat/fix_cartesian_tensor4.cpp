/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR4_CPP
#define CPPMAT_FIX_CARTESIAN_TENSOR4_CPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace cartesian {

// =================================================================================================
// constructors
// =================================================================================================

template<class X, size_t ND>
inline
tensor4<X,ND>::tensor4() : cppmat::tiny::array<X,4,ND,ND,ND,ND>()
{
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<class X, size_t ND>
template<typename U, typename V>
inline
tensor4<X,ND>::tensor4(const cppmat::tiny::array<U,4,ND,ND,ND,ND> &A) : cppmat::tiny::array<X,4,ND,ND,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t ND>
inline
tensor4<X,ND>::tensor4(const cppmat::cartesian::tensor4<X> &A) : cppmat::tiny::array<X,4,ND,ND,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t ND>
inline
tensor4<X,ND>::tensor4(const cppmat::view::cartesian::tensor4<X,ND> &A) : cppmat::tiny::array<X,4,ND,ND,ND,ND>(A)
{
}

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::I()
{
  tensor4<X,ND> out;

  out.setI();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::Irt()
{
  tensor4<X,ND> out;

  out.setIrt();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::Id()
{
  tensor4<X,ND> out;

  out.setId();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::Is()
{
  tensor4<X,ND> out;

  out.setIs();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::Isd()
{
  tensor4<X,ND> out;

  out.setIsd();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::II()
{
  tensor4<X,ND> out;

  out.setII();

  return out;
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
inline
size_t tensor4<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize: identity tensors
// =================================================================================================

template<class X, size_t ND>
inline
void tensor4<X,ND>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          if ( i == l and j == k )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void tensor4<X,ND>::setIrt()
{
  this->setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      for ( size_t k = 0 ; k < ND ; ++k )
        for ( size_t l = 0 ; l < ND ; ++l )
          if ( i == k and j == l )
            (*this)(i,j,k,l) = static_cast<X>(1);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void tensor4<X,ND>::setIs()
{
  (*this) = ( tensor4<X,ND>::I() + tensor4<X,ND>::Irt() ) / static_cast<X>(2);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void tensor4<X,ND>::setId()
{
  (*this) = tensor4<X,ND>::I() - tensor4<X,ND>::II()/static_cast<X>(ND);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void tensor4<X,ND>::setIsd()
{
  (*this) = tensor4<X,ND>::Is() - tensor4<X,ND>::II()/static_cast<X>(ND);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void tensor4<X,ND>::setII()
{
  this->setZero();

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

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::ddot(const tensor4<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor4<X,ND>::ddot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor4<X,ND>::ddot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor4<X,ND>::ddot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::T() const
{
  return cppmat::cartesian::T(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::RT() const
{
  return cppmat::cartesian::RT(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor4<X,ND>::LT() const
{
  return cppmat::cartesian::LT(*this);
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

