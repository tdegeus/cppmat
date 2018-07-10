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

// delete a specific item from a vector
template<class X> std::vector<X> del(const std::vector<X> &A, int    idx);
template<class X> std::vector<X> del(const std::vector<X> &A, size_t idx);

// return the indices that would sort the vector
template <typename X> std::vector<size_t> argsort(const std::vector<X> &v, bool ascending=true);

// convert vector items to string, and join these string together using the "join" string
template<class X> std::string to_string(const std::vector<X> &A, std::string join=", ");

// linearly spaced array
template <typename T = double> std::vector<T> linspace(T a, T b, size_t N);

// =================================================================================================

} // namespace ...

// =================================================================================================

// print operator
#ifndef CPPMAT_NOSTD
template<class X> std::ostream& operator<<(std::ostream& out, const std::vector<X>& src);
#endif

// =================================================================================================

#endif

