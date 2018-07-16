/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2S_HPP
#define CPPMAT_VAR_CARTESIAN_TENSOR2S_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// constructors
// =================================================================================================

template<typename X>
inline
tensor2s<X>::tensor2s(size_t nd) : cppmat::symmetric::matrix<X>(nd,nd)
{
  ND = nd;
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<typename X>
template<typename U, typename V>
inline
tensor2s<X>::tensor2s(const cppmat::symmetric::matrix<U> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<typename X>
inline
tensor2s<X>::tensor2s(const cppmat::diagonal::matrix<X> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<typename X>
template<size_t nd>
inline
tensor2s<X>::tensor2s(const cppmat::tiny::cartesian::tensor2s<X,nd> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<typename X>
template<size_t nd>
inline
tensor2s<X>::tensor2s(const cppmat::view::cartesian::tensor2s<X,nd> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<typename X>
inline
tensor2s<X> tensor2s<X>::Random(size_t nd, X lower, X upper)
{
  tensor2s<X> out(nd);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2s<X> tensor2s<X>::Arange(size_t nd)
{
  tensor2s<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2s<X> tensor2s<X>::Zero(size_t nd)
{
  tensor2s<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2s<X> tensor2s<X>::Ones(size_t nd)
{
  tensor2s<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2s<X> tensor2s<X>::Constant(size_t nd, X D)
{
  tensor2s<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::Copy(size_t nd, Iterator first)
{
  tensor2s<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  tensor2s<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::CopyDense(size_t nd, Iterator first)
{
  tensor2s<X> out(nd);

  out.setCopyDense(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<typename X>
template<typename Iterator>
inline
tensor2s<X> tensor2s<X>::CopyDense(size_t nd, Iterator first, Iterator last)
{
  tensor2s<X> out(nd);

  out.setCopyDense(first,last);

  return out;
}

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

template<typename X>
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

template<typename X>
inline
void tensor2s<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::symmetric::matrix<X>::resize(nd,nd);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
void tensor2s<X>::resize(size_t nd, const X &D)
{
  ND = nd;

  cppmat::symmetric::matrix<X>::resize(nd,nd, D);
}

// =================================================================================================
// dimensions
// =================================================================================================

template<typename X>
inline
size_t tensor2s<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// initialize: identity tensors
// =================================================================================================

template<typename X>
inline
void tensor2s<X>::setI()
{
  this->setZero();

  for ( size_t i = 0 ; i < ND ; ++i )
    (*this)(i,i) = static_cast<X>(1);
}

// =================================================================================================
// tensor products
// =================================================================================================

template<typename X>
inline
tensor2<X> tensor2s<X>::ddot(const tensor4<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
X tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<typename X>
inline
tensor2s<X> tensor2s<X>::T() const
{
  return cppmat::cartesian::T(*this);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<typename X>
inline
X tensor2s<X>::trace() const
{
  return cppmat::cartesian::trace(*this);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
X tensor2s<X>::det() const
{
  return cppmat::cartesian::det(*this);
}

// -------------------------------------------------------------------------------------------------

template<typename X>
inline
tensor2s<X> tensor2s<X>::inv() const
{
  return cppmat::cartesian::inv(*this);
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

