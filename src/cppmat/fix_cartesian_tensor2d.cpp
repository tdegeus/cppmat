/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR2D_CPP
#define CPPMAT_FIX_CARTESIAN_TENSOR2D_CPP

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
CPPMAT_INLINE
tensor2d<X,ND>::tensor2d() : cppmat::tiny::diagonal::matrix<X,ND,ND>()
{
}

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND>::tensor2d(const cppmat::tiny::diagonal::matrix<X,ND,ND> &A) : cppmat::tiny::diagonal::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND>::tensor2d(const cppmat::cartesian::tensor2d<X> &A) : cppmat::tiny::diagonal::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND>::tensor2d(const cppmat::view::cartesian::tensor2d<X,ND> &A) : cppmat::tiny::diagonal::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND> tensor2d<X,ND>::I()
{
  tensor2d<X,ND> out;

  out.setI();

  return out;
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
size_t tensor2d<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize: identity tensors
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
void tensor2d<X,ND>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    this->mData[i] = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2<X,ND> tensor2d<X,ND>::ddot(const tensor4<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
X tensor2d<X,ND>::ddot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
X tensor2d<X,ND>::ddot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
X tensor2d<X,ND>::ddot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor2<X,ND> tensor2d<X,ND>::dot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor2<X,ND> tensor2d<X,ND>::dot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND> tensor2d<X,ND>::dot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
vector<X,ND> tensor2d<X,ND>::dot(const vector<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor4<X,ND> tensor2d<X,ND>::dyadic(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor4<X,ND> tensor2d<X,ND>::dyadic(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor4<X,ND> tensor2d<X,ND>::dyadic(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND> tensor2d<X,ND>::T() const
{
  return cppmat::cartesian::T(*this);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X, size_t ND>
CPPMAT_INLINE
X tensor2d<X,ND>::trace() const
{
  return cppmat::cartesian::trace(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
X tensor2d<X,ND>::det() const
{
  return cppmat::cartesian::det(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
CPPMAT_INLINE
tensor2d<X,ND> tensor2d<X,ND>::inv() const
{
  return cppmat::cartesian::inv(*this);
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

