/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_VECTOR_CPP
#define CPPMAT_VAR_CARTESIAN_VECTOR_CPP

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
vector<X>::vector(size_t nd) : cppmat::vector<X>(nd)
{
  ND = nd;
}

// =================================================================================================
// constructors: copy from parent
// =================================================================================================

template<class X>
inline
vector<X>::vector(const cppmat::array<X> &A) : cppmat::vector<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X>
inline
vector<X>::vector(const std::vector<X> &A) : cppmat::vector<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from fixed size
// =================================================================================================

template<class X>
template<size_t nd>
inline
vector<X>::vector(const cppmat::tiny::cartesian::vector<X,nd> &A) : cppmat::array<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X>
template<size_t nd>
inline
vector<X>::vector(const cppmat::view::cartesian::vector<X,nd> &A) : cppmat::array<X>(A)
{
  ND = this->mShape[0];
}

// =================================================================================================
// named constructors
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::Random(size_t nd, X lower, X upper)
{
  vector<X> out(nd);

  out.setRandom(lower, upper);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Arange(size_t nd)
{
  vector<X> out(nd);

  out.setArange();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Zero(size_t nd)
{
  vector<X> out(nd);

  out.setZero();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Ones(size_t nd)
{
  vector<X> out(nd);

  out.setOnes();

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Constant(size_t nd, X D)
{
  vector<X> out(nd);

  out.setConstant(D);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Copy(const std::vector<X> &D)
{
  vector<X> out(D.size());

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::Copy(size_t nd, const std::vector<X> &D)
{
  vector<X> out(nd);

  out.setCopy(D.begin(), D.end());

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(size_t nd, Iterator first)
{
  vector<X> out(nd);

  out.setCopy(first);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(size_t nd, Iterator first, Iterator last)
{
  vector<X> out(nd);

  out.setCopy(first,last);

  return out;
}

// -------------------------------------------------------------------------------------------------

template<class X>
template<typename Iterator>
inline
vector<X> vector<X>::Copy(Iterator first, Iterator last)
{
  vector<X> out(last-first);

  out.setCopy(first,last);

  return out;
}

// =================================================================================================
// resize
// =================================================================================================

template<class X>
inline
void vector<X>::resize(size_t nd)
{
  ND = nd;

  cppmat::vector<X>::resize(nd);
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X>
inline
size_t vector<X>::ndim() const
{
  return ND;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X>
inline
vector<X> vector<X>::dot(const tensor2<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X vector<X>::dot(const vector<X> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::cross(const vector<X> &B) const
{
  return cppmat::cartesian::cross(*this, B);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X vector<X>::length() const
{
  return cppmat::cartesian::length(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void vector<X>::setUnitLength()
{
  X C = this->length();

  if ( C <= static_cast<X>(0) ) return;

  for ( size_t i = 0 ; i < this->mSize ; ++i )
    this->mData[i] /= C;
}

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

