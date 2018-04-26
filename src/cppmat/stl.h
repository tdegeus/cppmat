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

template<class X> std::vector<X> sort_pmodulo(const std::vector<X> &A, X n, bool reverse=true);

// -------------------------------------------------------------------------------------------------

template<class X> std::string to_string(const std::vector<X> &A);

// -------------------------------------------------------------------------------------------------

template<class X> X                        abs(X                               A);
template<class X> matrix<X>                abs(const matrix<X>                &A);
template<class X> cartesian2d::tensor2 <X> abs(const cartesian2d::tensor2 <X> &A);
template<class X> cartesian2d::tensor2s<X> abs(const cartesian2d::tensor2s<X> &A);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

