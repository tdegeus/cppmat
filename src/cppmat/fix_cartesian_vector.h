/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_VECTOR_H
#define CPPMAT_FIX_CARTESIAN_VECTOR_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace cartesian {

// =================================================================================================
// cppmat::tiny::cartesian::vector
// =================================================================================================

template<class X, size_t ND>
class vector : public cppmat::tiny::vector<X,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  vector();

  // constructor: copy from parent
  vector(const cppmat::tiny::array<X,1,ND> &A);

  // constructor: copy from other classes
  vector(const std::vector<X> &A);

  // constructor: copy from dynamic size
  vector(const cppmat::cartesian::vector<X> &A);

  // constructor: copy from view
  vector(const cppmat::view::cartesian::vector<X,ND> &A);

  // get dimensions
  size_t ndim() const;

  // tensor products / operations
  X             dot   (const vector  <X,ND> &B) const; // dot    product: C   = A_i*B_i
  vector <X,ND> dot   (const tensor2 <X,ND> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X,ND> dot   (const tensor2s<X,ND> &B) const; // dot    product: C_j = A_i*B_ij
  vector <X,ND> dot   (const tensor2d<X,ND> &B) const; // dot    product: C_j = A_i*B_ij
  tensor2<X,ND> dyadic(const vector  <X,ND> &B) const; // dyadic product: C_ij = A_i*B_j
  vector <X,ND> cross (const vector  <X,ND> &B) const; // cross  product (only in 3D)
  X             length()                        const; // sqrt(sum(pow(A_i,2.)))
  void          setUnitLength();                       // A_i /= A.length()

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

