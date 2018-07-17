/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_REGULAR_VECTOR_H
#define CPPMAT_FIX_REGULAR_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

// =================================================================================================
// cppmat::tiny::vector
// =================================================================================================

template<typename X, size_t N>
class vector : public cppmat::tiny::array<X,1,N>
{
private:
  size_t mIstore=0;

public:

  // constructor: allocate, don't initialize
  vector();

  // constructor: copy from parent (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  vector(const cppmat::tiny::array<U,1,N> &A);

  // constructor: copy from other class
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  vector(const std::vector<U> &A);

  // constructor: copy from {...}
  vector(const std::initializer_list<X> &A);

  // constructor: copy from dynamic size
  vector(const cppmat::vector<X> &A);

  // constructor: copy from view
  vector(const cppmat::view::vector<X,N> &A);

  // forward difference (x0, x1-x0, x2-x1, ...)
  vector<X,N> diff() const;

  // STL-like access
  void push_back(const X &value);
  void clear();


};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

