/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_VECTOR_CPP
#define CPPMAT_FIX_CARTESIAN_VECTOR_CPP

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
vector<X,ND>::vector() : cppmat::tiny::vector<X,ND>()
{
}

// =================================================================================================
// constructors: copy from parent (with different type)
// =================================================================================================

template<class X, size_t ND>
template<typename U, typename V>
inline
vector<X,ND>::vector(const cppmat::tiny::array<U,1,ND> &A) : cppmat::tiny::vector<X,ND>(A)
{
}

// =================================================================================================
// constructors: copy from other class
// =================================================================================================

template<class X, size_t ND>
inline
vector<X,ND>::vector(const std::vector<X> &A) : cppmat::tiny::vector<X,ND>(A)
{
}

// =================================================================================================
// constructors: copy from dynamic size
// =================================================================================================

template<class X, size_t ND>
inline
vector<X,ND>::vector(const cppmat::cartesian::vector<X> &A) : cppmat::tiny::vector<X,ND>(A)
{
}

// =================================================================================================
// constructors: copy from view
// =================================================================================================

template<class X, size_t ND>
inline
vector<X,ND>::vector(const cppmat::view::cartesian::vector<X,ND> &A) : cppmat::tiny::vector<X,ND>(A)
{
}

// =================================================================================================
// dimensions
// =================================================================================================

template<class X, size_t ND>
inline
size_t vector<X,ND>::ndim() const
{
  return ND;
}

// =================================================================================================
// tensor products
// =================================================================================================

template<class X, size_t ND>
inline
vector<X,ND> vector<X,ND>::dot(const tensor2<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
vector<X,ND> vector<X,ND>::dot(const tensor2s<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
vector<X,ND> vector<X,ND>::dot(const tensor2d<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
X vector<X,ND>::dot(const vector<X,ND> &B) const
{
  return cppmat::cartesian::dot(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
tensor2<X,ND> vector<X,ND>::dyadic(const vector<X,ND> &B) const
{
  return cppmat::cartesian::dyadic(*this, B);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
vector<X,ND> vector<X,ND>::cross(const vector<X,ND> &B) const
{
  return cppmat::cartesian::cross(*this, B);
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X, size_t ND>
inline
X vector<X,ND>::length() const
{
  return cppmat::cartesian::length(*this);
}

// -------------------------------------------------------------------------------------------------

template<class X, size_t ND>
inline
void vector<X,ND>::setUnitLength()
{
  X C = this->length();

  if ( C <= static_cast<X>(0) ) return;

  for ( size_t i = 0 ; i < this->mSize ; ++i )
    this->mData[i] /= C;
}

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

