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

template<class X> std::string to_string(const std::vector<X> &A);

// -------------------------------------------------------------------------------------------------

template<class X> array<X>  abs(const array<X>  &A);
template<class X> matrix<X> abs(const matrix<X> &A);
template<class X> vector<X> abs(const vector<X> &A);

template<class X> periodic::array<X>  abs(const periodic::array<X>  &A);
template<class X> periodic::matrix<X> abs(const periodic::matrix<X> &A);
template<class X> periodic::vector<X> abs(const periodic::vector<X> &A);

template<class X, size_t m, size_t n> tiny::matrix<X,m,n> abs(const tiny::matrix<X,m,n> &A);
template<class X, size_t n>           tiny::vector<X,n>   abs(const tiny::vector<X,n>   &A);

template<class X> cartesian2d::tensor2 <X> abs(const cartesian2d::tensor2 <X> &A);
template<class X> cartesian2d::tensor2s<X> abs(const cartesian2d::tensor2s<X> &A);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif

