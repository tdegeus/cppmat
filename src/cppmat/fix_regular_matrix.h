/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_MATRIX_H
#define CPPMAT_FIX_REGULAR_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::matrix
// =================================================================================================

template<class X, size_t M, size_t N>
class matrix : public cppmat::tiny::array<X,2,M,N>
{
public:

  // avoid name-hiding (also use the already defined overloads)
  using cppmat::tiny::array<X,2,M,N>::operator*=;
  using cppmat::tiny::array<X,2,M,N>::operator/=;
  using cppmat::tiny::array<X,2,M,N>::operator+=;
  using cppmat::tiny::array<X,2,M,N>::operator-=;

public:

  // constructor: allocate, don't initialize
  matrix();

  // constructor: copy
  matrix(const cppmat::tiny::           array <X,2,M,N> &A);
  matrix(const cppmat::tiny::symmetric::matrix<X,  M,N> &A);
  matrix(const cppmat::tiny::diagonal ::matrix<X,  M,N> &A);


  // extra arithmetic operators
  matrix<X,M,N>& operator*= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);
  matrix<X,M,N>& operator/= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);
  matrix<X,M,N>& operator+= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);
  matrix<X,M,N>& operator-= (const cppmat::tiny::symmetric::matrix<X,M,N> &B);
  matrix<X,M,N>& operator*= (const cppmat::tiny::diagonal ::matrix<X,M,N> &B);
  matrix<X,M,N>& operator+= (const cppmat::tiny::diagonal ::matrix<X,M,N> &B);
  matrix<X,M,N>& operator-= (const cppmat::tiny::diagonal ::matrix<X,M,N> &B);
};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

