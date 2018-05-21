/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR2_H
#define CPPMAT_VAR_CARTESIAN_TENSOR2_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// cppmat::cartesian::tensor2
// =================================================================================================

template<class X>
class tensor2 : public cppmat::matrix<X>
{
protected:

  size_t ND=0;

public:

  // constructor
  tensor2() = default;

  // constructor: allocate, don't initialize
  tensor2(size_t nd);

  // constructor: copy
  tensor2(const cppmat::           array <X> &A);
  tensor2(const cppmat::           matrix<X> &A);
  tensor2(const cppmat::symmetric::matrix<X> &A);
  tensor2(const cppmat::diagonal ::matrix<X> &A);

  // named constructor: initialize
  static tensor2<X> Random  (size_t nd, X lower=(X)0, X upper=(X)1);
  static tensor2<X> Arange  (size_t nd);
  static tensor2<X> Zero    (size_t nd);
  static tensor2<X> Ones    (size_t nd);
  static tensor2<X> Constant(size_t nd, X D);
  static tensor2<X> I       (size_t nd);

  // named constructor: copy
  static tensor2<X> Copy(size_t nd, const std::vector<X> &D);

  // named constructor: copy
  template<typename Iterator> static tensor2<X> Copy(size_t nd, Iterator first);
  template<typename Iterator> static tensor2<X> Copy(size_t nd, Iterator first, Iterator last);

  // resize
  void resize(size_t nd);

  // number of dimensions (== shape[0]...)
  size_t ndim() const;

  // initialize
  void setI();

  // tensor products / operations
  tensor2 <X> dot   (const tensor2 <X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> dot   (const tensor2s<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  tensor2 <X> dot   (const tensor2d<X> &B) const; // single contract.: C_ik   = A_ij * B_jk
  vector  <X> dot   (const vector  <X> &B) const; // single contract.: C_i    = A_ij * B_j
  tensor2 <X> ddot  (const tensor4 <X> &B) const; // double contract.: C_kl   = A_ij * B_jikl
  X           ddot  (const tensor2 <X> &B) const; // double contract.: C      = A_ij * B_ji
  X           ddot  (const tensor2s<X> &B) const; // double contract.: C      = A_ij * B_ji
  X           ddot  (const tensor2d<X> &B) const; // double contract.: C      = A_ij * B_ji
  tensor4 <X> dyadic(const tensor2 <X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> dyadic(const tensor2s<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor4 <X> dyadic(const tensor2d<X> &B) const; // dyadic product  : C_ijkl = A_ij * B_kl
  tensor2 <X> T     ()                     const; // transpose       : C_ij   = A_ji
  X           trace ()                     const; // trace           : A_ii
  X           det   ()                     const; // determinant (only in 2D/3D)
  tensor2 <X> inv   ()                     const; // inverse     (only in 2D/3D)

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

