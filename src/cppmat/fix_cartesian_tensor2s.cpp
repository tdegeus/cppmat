/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR2S_CPP
#define CPPMAT_FIX_CARTESIAN_TENSOR2S_CPP

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
tensor2s<X,ND>::tensor2s() : cppmat::tiny::symmetric::matrix<X,ND,ND>()
{
}

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s(const cppmat::tiny::symmetric::matrix<X,ND,ND> &A) : cppmat::tiny::symmetric::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s(const cppmat::tiny::diagonal::matrix<X,ND,ND> &A) : cppmat::tiny::symmetric::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s(const cppmat::cartesian::tensor2s<X> &A) : cppmat::tiny::symmetric::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND>::tensor2s(const cppmat::view::cartesian::tensor2s<X,ND> &A) : cppmat::tiny::symmetric::matrix<X,ND,ND>(A)
{
}

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND> tensor2s<X,ND>::I()
{
  tensor2s<X,ND> out;

  out.setI();

  return out;
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
inline
size_t tensor2s<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize: identity tensors
// =================================================================================================

template<class X, size_t ND>
inline
void tensor2s<X,ND>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X, size_t ND>
inline
tensor2<X,ND> tensor2s<X,ND>::ddot(const tensor4<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X tensor2s<X,ND>::ddot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X tensor2s<X,ND>::ddot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X tensor2s<X,ND>::ddot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor2s<X,ND>::dot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor2s<X,ND>::dot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> tensor2s<X,ND>::dot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
vector<X,ND> tensor2s<X,ND>::dot(const vector<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor2s<X,ND>::dyadic(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor2s<X,ND>::dyadic(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor4<X,ND> tensor2s<X,ND>::dyadic(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X, size_t ND>
inline
tensor2s<X,ND> tensor2s<X,ND>::T() const
{
  return cppmat::cartesian::T(*this);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X, size_t ND>
inline
X tensor2s<X,ND>::trace() const
{
  return cppmat::cartesian::trace(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X tensor2s<X,ND>::det() const
{
  return cppmat::cartesian::det(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2s<X,ND> tensor2s<X,ND>::inv() const
{
  return cppmat::cartesian::inv(*this);
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

