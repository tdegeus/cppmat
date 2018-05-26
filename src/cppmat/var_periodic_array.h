/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_PERIODIC_ARRAY_H
#define CPPMAT_VAR_PERIODIC_ARRAY_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

// =================================================================================================
// cppmat::periodic::array
// =================================================================================================

template<class X>
class array : public cppmat::array<X>
{
protected:

  int mShapeI  [cppmat::array<X>::MAX_DIM];  // == mShape   with int's
  int mStridesI[cppmat::array<X>::MAX_DIM];  // == mStrides with int's

public:

  // avoid name-hiding (also use the already defined overloads)
  using cppmat::array<X>::operator();
  using cppmat::array<X>::compress;
  using cppmat::array<X>::item;

public:

  // constructor: default
  array() = default;

  // constructor: allocate, don't initialize
  array(const std::vector<size_t> &shape);

  // constructor: copy from parent
  array(const cppmat::array<X> &A);

  // constructor: copy from own class
  array(const cppmat::periodic::array<X> &A);

  // constructor: copy from fixed size
  template<size_t rank, size_t i, size_t j, size_t k, size_t l, size_t m, size_t n>
  array(const cppmat::tiny::periodic::array<X,rank,i,j,k,l,m,n> &A);

  // constructor: copy from view
  template<size_t rank, size_t i, size_t j, size_t k, size_t l, size_t m, size_t n>
  array(const cppmat::view::periodic::array<X,rank,i,j,k,l,m,n> &A);

  // named constructor: initialize
  static array<X> Random  (const std::vector<size_t> &shape, X lower=(X)0, X upper=(X)1);
  static array<X> Arange  (const std::vector<size_t> &shape);
  static array<X> Zero    (const std::vector<size_t> &shape);
  static array<X> Ones    (const std::vector<size_t> &shape);
  static array<X> Constant(const std::vector<size_t> &shape, X D);
  static array<X> Copy    (const std::vector<size_t> &shape, const std::vector<X> &D);

  // named constructor: copy
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first);
  template<typename It> static array<X> Copy(const std::vector<size_t> &shape, It first, It last);

  // resize
  void resize (const std::vector<size_t> &shape);
  void reshape(const std::vector<size_t> &shape);

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

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

