/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR2D_H
#define CPPMAT_FIX_CARTESIAN_TENSOR2D_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace cartesian {

// =================================================================================================
// cppmat::tiny::cartesian::tensor2d
// =================================================================================================

template<class X, size_t ND>
class tensor2d : public cppmat::tiny::diagonal::matrix<X,ND,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  tensor2d();

  // constructor: copy from parent
  tensor2d(const cppmat::tiny::diagonal::matrix<X,ND,ND> &A);

  // constructor: copy from dynamic size
  tensor2d(const cppmat::cartesian::tensor2d<X> &A);

  // constructor: copy from view
  tensor2d(const cppmat::view::cartesian::tensor2d<X,ND> &A);

  // named constructor: initialize
  static tensor2d<X,ND> I();

  // get dimensions
  size_t ndim() const;

  // initialize
  void setI();

  // tensor products / operations
  tensor2 <X,ND> dot   (const tensor2 <X,ND> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X,ND> dot   (const tensor2s<X,ND> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2d<X,ND> dot   (const tensor2d<X,ND> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X,ND> dot   (const vector  <X,ND> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X,ND> ddot  (const tensor4 <X,ND> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X              ddot  (const tensor2 <X,ND> &B) const; // double contract.: C      = A_ij * B_ji
  X              ddot  (const tensor2s<X,ND> &B) const; // double contract.: C      = A_ij * B_ji
  X              ddot  (const tensor2d<X,ND> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X,ND> dyadic(const tensor2 <X,ND> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X,ND> dyadic(const tensor2s<X,ND> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X,ND> dyadic(const tensor2d<X,ND> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2d<X,ND> T     ()                        const; // transpose       : C_ij   = A_ji
  X              trace ()                        const; // trace           : A_ii
  X              det   ()                        const; // determinant (only in 2D/3D)
  tensor2d<X,ND> inv   ()                        const; // inverse     (only in 2D/3D)

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

