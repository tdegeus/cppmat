/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_CARTESIAN_TENSOR4_H
#define CPPMAT_VAR_CARTESIAN_TENSOR4_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian {

// =================================================================================================
// cppmat::cartesian::tensor4
// =================================================================================================

template<typename X>
class tensor4 : public cppmat::array<X>
{
protected:

  // local variables
  size_t ND=0; // number of dimensions (== mShape[0] == mShape[1] == ...)

private:

  // hide functions
  using cppmat::array<X>::chrank;

public:

  // constructor: default
  tensor4() = default;

  // constructor: allocate, don't initialize
  tensor4(size_t nd);

  // constructor: copy from parent (with different type)
  template<typename U, typename=typename std::enable_if<std::is_convertible<U,X>::value>::type>
  tensor4(const cppmat::array<U> &A);

  // constructor: copy from fixed size
  template<size_t nd> tensor4(const cppmat::tiny::cartesian::tensor4<X,nd> &A);

  // constructor: copy from view
  template<size_t nd> tensor4(const cppmat::view::cartesian::tensor4<X,nd> &A);

  // named constructor: initialize
  static tensor4<X> Random  (size_t nd, X lower=(X)0, X upper=(X)1);
  static tensor4<X> Arange  (size_t nd);
  static tensor4<X> Zero    (size_t nd);
  static tensor4<X> Ones    (size_t nd);
  static tensor4<X> Constant(size_t nd, X D);
  static tensor4<X> Copy    (size_t nd, const std::vector<X> &D);
  static tensor4<X> I       (size_t nd);
  static tensor4<X> Irt     (size_t nd);
  static tensor4<X> Id      (size_t nd);
  static tensor4<X> Is      (size_t nd);
  static tensor4<X> Isd     (size_t nd);
  static tensor4<X> II      (size_t nd);

  // named constructor: copy
  template<typename Iterator> static tensor4<X> Copy(size_t nd, Iterator first);
  template<typename Iterator> static tensor4<X> Copy(size_t nd, Iterator first, Iterator last);

  // resize
  void resize(size_t nd);
  void resize(size_t nd, const X &D);

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
  tensor4<X> ddot(const tensor4 <X> &B) const; // double contract.: C_ijmn = A_ijkl * B_lkmn
  tensor2<X> ddot(const tensor2 <X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> ddot(const tensor2s<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor2<X> ddot(const tensor2d<X> &B) const; // double contract.: C_ij   = A_ijkl * B_lk
  tensor4<X> T   ()                     const; // transposition   : C_lkji = A_ijkl
  tensor4<X> RT  ()                     const; // transposition   : C_ijlk = A_ijkl
  tensor4<X> LT  ()                     const; // transposition   : C_jikl = A_ijkl

};

// =================================================================================================

}} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

