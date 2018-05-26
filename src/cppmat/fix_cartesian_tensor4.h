/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_FIX_CARTESIAN_TENSOR4_H
#define CPPMAT_FIX_CARTESIAN_TENSOR4_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace cartesian {

// =================================================================================================
// cppmat::tiny::cartesian::tensor4
// =================================================================================================

template<class X, size_t ND>
class tensor4 : public cppmat::tiny::array<X,4,ND,ND,ND,ND>
{
  static_assert( ND > 0, "Number of dimensions must positive" );

public:

  // constructor: allocate, don't initialize
  tensor4();

  // constructor: copy from parent
  tensor4(const cppmat::tiny::array<X,4,ND,ND,ND,ND> &A);

  // constructor: copy from dynamic size
  tensor4(const cppmat::cartesian::tensor4<X> &A);

  // constructor: copy from view
  tensor4(const cppmat::view::cartesian::tensor4<X,ND> &A);

  // named constructor: initialize
  static tensor4<X,ND> I  ();
  static tensor4<X,ND> Irt();
  static tensor4<X,ND> Id ();
  static tensor4<X,ND> Is ();
  static tensor4<X,ND> Isd();
  static tensor4<X,ND> II ();

  // get dimensions
  size_t ndim() const;

  // initialize
  void setI();
  void setIrt();
  void setId();
  void setIs();
  void setIsd();
  void setII();

  // tensor products / operations
  tensor4<X,ND> ddot(const tensor4 <X,ND> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2<X,ND> ddot(const tensor2 <X,ND> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X,ND> ddot(const tensor2s<X,ND> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X,ND> ddot(const tensor2d<X,ND> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor4<X,ND> T   ()                        const; // transposition   : C_lkji = A_ijkl
  tensor4<X,ND> RT  ()                        const; // transposition   : C_ijlk = A_ijkl
  tensor4<X,ND> LT  ()                        const; // transposition   : C_jikl = A_ijkl

};

// =================================================================================================

}}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

