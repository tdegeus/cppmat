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

template<class X>
inline
tensor2s<X>::tensor2s(const cppmat::diagonal::matrix<X> &A) : cppmat::symmetric::matrix<X>(A)
{
  ND = this->N;
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
tensor2s<X> tensor2s<X>::Random(size_t nd, X lower, X upper)
{
  tensor2s<X> out(nd);

  out.setRandom(lower, upper);

  return out;
}

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

// =================================================================================================
// named constructors: identity tensors
// =================================================================================================

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
// initialize: identity tensors
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
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::ddot(const tensor2<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::ddot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::ddot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::ddot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> tensor2s<X>::dot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> tensor2s<X>::dot(const vector<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor4<X> tensor2s<X>::dyadic(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// =================================================================================================
// transpositions
// =================================================================================================

template<class X>
inline
tensor2s<X> tensor2s<X>::T() const
{
  return cppmat::cartesian::T(*this);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X tensor2s<X>::trace() const
{
  return cppmat::cartesian::trace(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X tensor2s<X>::det() const
{
  return cppmat::cartesian::det(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2s<X> tensor2s<X>::inv() const
{
  return cppmat::cartesian::inv(*this);
}

// =================================================================================================

}} // namespace ...

#endif

