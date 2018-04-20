/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_STL_H
#define CPPMAT_STL_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================

template<class X> std::vector<X> del(const std::vector<X> &A, int    idx);
template<class X> std::vector<X> del(const std::vector<X> &A, size_t idx);

// -------------------------------------------------------------------------------------------------

template<class X> X                        abs(X                               A);
template<class X> matrix<X>                abs(const matrix<X>                &A);
template<class X> cartesian2d::tensor2 <X> abs(const cartesian2d::tensor2 <X> &A);
template<class X> cartesian2d::tensor2s<X> abs(const cartesian2d::tensor2s<X> &A);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

