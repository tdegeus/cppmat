/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2D_H
#define CPPMAT_VAR_CARTESIAN_TENSOR2D_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// cppmat::cartesian::tensor2d
// =================================================================================================

template<class X>
class tensor2d : public cppmat::diagonal::matrix<X>
{
protected:

  // local variables
  size_t ND=0; // number of dimensions (== mShape[0] == mShape[1])

public:

  // constructor: default
  tensor2d();

  // constructor: allocate, don't initialize
  tensor2d(size_t nd);

  // constructor: copy
  tensor2d(const cppmat::diagonal::matrix<X> &A);

  // named constructor: initialize
  static tensor2d<X> Random  (size_t nd, X lower=(X)0, X upper=(X)1);
  static tensor2d<X> Arange  (size_t nd);
  static tensor2d<X> Zero    (size_t nd);
  static tensor2d<X> Ones    (size_t nd);
  static tensor2d<X> Constant(size_t nd, X D);
  static tensor2d<X> I       (size_t nd);

  // named constructor: copy
  template<typename Iterator> static tensor2d<X> Copy     (size_t nd, Iterator first);
  template<typename Iterator> static tensor2d<X> Copy     (size_t nd, Iterator first, Iterator last);
  template<typename Iterator> static tensor2d<X> CopyDense(size_t nd, Iterator first);
  template<typename Iterator> static tensor2d<X> CopyDense(size_t nd, Iterator first, Iterator last);

  // resize
  void resize(size_t nd);

  // get dimensions
  size_t ndim() const;

  // initialize
  void setI();

  // tensor products / operations
  tensor2 <X> dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2d<X> dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X> dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X> ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X           ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X           ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X           ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X> dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2d<X> T     ()                     const; // transpose       : C_ij   = A_ji
  X           trace ()                     const; // trace           : A_ii
  X           det   ()                     const; // determinant (only in 2D/3D)
  tensor2d<X> inv   ()                     const; // inverse     (only in 2D/3D)

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

