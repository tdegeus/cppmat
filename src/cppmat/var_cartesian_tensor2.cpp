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

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X>
inline
tensor2<X>::tensor2(const cppmat::array<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X>
inline
tensor2<X>::tensor2(const cppmat::symmetric::matrix<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->mShape[0];
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X>::tensor2(const cppmat::diagonal::matrix<X> &A) : cppmat::matrix<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<class X>
template<size_t nd>
inline
tensor2<X>::tensor2(const cppmat::tiny::cartesian::tensor2<X,nd> &A) : cppmat::matrix<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X>
template<size_t nd>
inline
tensor2<X>::tensor2(const cppmat::view::cartesian::tensor2<X,nd> &A) : cppmat::matrix<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2<X>::Random(size_t nd, X lower, X upper)
{
  tensor2<X> out(nd);

  out.setRandom(lower, upper);

  return out;
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
inline
tensor2<X> tensor2<X>::Copy(size_t nd, const std::vector<X> &D)
{
  tensor2<X> out(nd);

  out.setCopy(D.begin(), D.end());

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

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

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
// initialize: identity tensors
// =================================================================================================

template<class X>
inline
void tensor2<X>::setI()
{
  this->setZero();

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
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::ddot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::dot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> tensor2<X>::dot(const vector<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2<X>::dyadic(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
tensor2<X> tensor2<X>::T() const
{
  return cppmat::cartesian::T(*this);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X tensor2<X>::trace() const
{
  return cppmat::cartesian::trace(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2<X>::det() const
{
  return cppmat::cartesian::det(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2<X>::inv() const
{
  return cppmat::cartesian::inv(*this);
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

