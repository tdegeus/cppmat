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

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const cppmat::array<X> &A) : cppmat::vector<X>(A)
{
  ND = this->N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const cppmat::vector<X> &A) : cppmat::vector<X>(A)
{
  ND = this->N;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X>::vector(const std::vector<X> &D) : cppmat::vector<X>(D)
{
  ND = this->N;
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
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::dot(const tensor2s<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(j) += (*this)(i) * B(i,j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::dot(const tensor2d<X> &B) const
{
  assert( ndim() == B.ndim() );

  vector<X> C = vector<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    C(i) += (*this)(i) * B(i,i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
X vector<X>::dot(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < ND ; ++i )
    C += (*this)(i) * B(i);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
tensor2<X> vector<X>::dyadic(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  tensor2<X> C = tensor2<X>::Zero(ND);

  for ( size_t i = 0 ; i < ND ; ++i )
    for ( size_t j = 0 ; j < ND ; ++j )
      C(i,j) += (*this)(i) * B(j);

  return C;
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
vector<X> vector<X>::cross(const vector<X> &B) const
{
  assert( ndim() == B.ndim() );

  if ( ND != 3 )
    throw std::runtime_error("'cross' only implemented in 3D");

  vector<X> C(3);

  C[0] =                     this->mData[1]*B[2]-B[1]*this->mData[2] ;
  C[1] = static_cast<X>(-1)*(this->mData[0]*B[2]-B[0]*this->mData[2]);
  C[2] =                     this->mData[0]*B[1]-B[0]*this->mData[1] ;

  return C;
}

// =================================================================================================
// miscellaneous tensor operations
// =================================================================================================

template<class X>
inline
X vector<X>::length() const
{
  X C = static_cast<X>(0);

  for ( size_t i = 0 ; i < this->mSize ; ++i )
    C += std::pow(this->mData[i],2.);

  return std::sqrt(C);
}

// -------------------------------------------------------------------------------------------------

template<class X>
inline
void vector<X>::setUnitLength()
{
  X C = length();

  if ( C <= static_cast<X>(0) ) return;

  for ( size_t i = 0 ; i < this->mSize ; ++i )
    this->mData[i] /= C;
}

// =================================================================================================

}} // namespace ...

#endif

