/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_PERIODIC_ARRAY_H
#define CPPMAT_FIX_PERIODIC_ARRAY_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace periodic {

// =================================================================================================
// cppmat::tiny::periodic::array
// =================================================================================================

template<class X, size_t RANK, size_t I, size_t J, size_t K, size_t L, size_t M, size_t N>
class array : public cppmat::tiny::array<X,RANK,I,J,K,L,M,N>
{
protected:

  int mShapeI  [cppmat::tiny::array<X,RANK,I,J,K,L,M,N>::MAX_DIM];  // == mShape   with int's
  int mStridesI[cppmat::tiny::array<X,RANK,I,J,K,L,M,N>::MAX_DIM];  // == mStrides with int's

public:

  // avoid name-hiding (also use the already defined overloads)
  using cppmat::tiny::array<X,RANK,I,J,K,L,M,N>::operator();
  using cppmat::tiny::array<X,RANK,I,J,K,L,M,N>::compress;
  using cppmat::tiny::array<X,RANK,I,J,K,L,M,N>::item;

public:

  // constructor: allocate, don't initialize
  array();

  // constructor: copy from parent
  array(const cppmat::tiny::array<X,RANK,I,J,K,L,M,N> &A);

  // constructor: copy from own class
  array(const cppmat::tiny::periodic::array<X,RANK,I,J,K,L,M,N> &A);

  // constructor: copy from dynamic size
  array(const cppmat::periodic::array<X> &A);

  // constructor: copy from view
  array(const cppmat::view::periodic::array<X,RANK,I,J,K,L,M,N> &A);

  // index operators: access using array-indices
  X&       operator()(int a);
  const X& operator()(int a) const;
  X&       operator()(int a, int b);
  const X& operator()(int a, int b) const;
  X&       operator()(int a, int b, int c);
  const X& operator()(int a, int b, int c) const;
  X&       operator()(int a, int b, int c, int d);
  const X& operator()(int a, int b, int c, int d) const;
  X&       operator()(int a, int b, int c, int d, int e);
  const X& operator()(int a, int b, int c, int d, int e) const;
  X&       operator()(int a, int b, int c, int d, int e, int f);
  const X& operator()(int a, int b, int c, int d, int e, int f) const;

  // index operators: access using iterator
  // N.B. the iterator points to list of array-indices (a,b,c,...)
  template<class Iterator> X&       at(Iterator first, Iterator last);
  template<class Iterator> const X& at(Iterator first, Iterator last) const;

  // index operators: array-indices -> plain storage (a,b,c,... -> i)
  size_t compress(int a) const;
  size_t compress(int a, int b) const;
  size_t compress(int a, int b, int c) const;
  size_t compress(int a, int b, int c, int d) const;
  size_t compress(int a, int b, int c, int d, int e) const;
  size_t compress(int a, int b, int c, int d, int e, int f) const;

  // iterator to specific entry: access using array-indices
  auto item(int a);
  auto item(int a) const;
  auto item(int a, int b);
  auto item(int a, int b) const;
  auto item(int a, int b, int c);
  auto item(int a, int b, int c) const;
  auto item(int a, int b, int c, int d);
  auto item(int a, int b, int c, int d) const;
  auto item(int a, int b, int c, int d, int e);
  auto item(int a, int b, int c, int d, int e) const;
  auto item(int a, int b, int c, int d, int e, int f);
  auto item(int a, int b, int c, int d, int e, int f) const;

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

