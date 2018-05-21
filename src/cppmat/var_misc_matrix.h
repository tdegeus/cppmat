/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_VAR_MISC_MATRIX_H
#define CPPMAT_VAR_MISC_MATRIX_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// =================================================================================================
// extra external arithmetic operators -> cppmat::matrix
// =================================================================================================

template<class X>
inline
cppmat::matrix<X> operator* (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator/ (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator* (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator/ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
);

// =================================================================================================
// extra external arithmetic operators -> cppmat::symmetric::matrix
// =================================================================================================

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const X &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const X &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator+ (
  const X &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::symmetric::matrix<X> operator- (
  const X &A,
  const cppmat::diagonal::matrix<X> &B
);

// =================================================================================================
// extra external arithmetic operators -> cppmat::diagonal::matrix
// =================================================================================================

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::diagonal::matrix<X> operator/ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::matrix<X> &B
);

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::diagonal::matrix<X> operator/ (
  const cppmat::diagonal::matrix<X> &A,
  const cppmat::symmetric::matrix<X> &B
);

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

template<class X>
inline
cppmat::diagonal::matrix<X> operator* (
  const cppmat::symmetric::matrix<X> &A,
  const cppmat::diagonal::matrix<X> &B
);

// =================================================================================================

#endif

