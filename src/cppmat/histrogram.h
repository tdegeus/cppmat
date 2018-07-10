/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_HISTOGRAM_H
#define CPPMAT_HISTOGRAM_H

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================

template<typename X>
std::tuple<std::vector<double>, std::vector<double>> histogram(
  const std::vector<X> &data, size_t bins, bool density=false, bool return_edges=false);

// -------------------------------------------------------------------------------------------------

template<typename X>
std::tuple<std::vector<double>, std::vector<double>> histogram_uniform(
  const std::vector<X> &data, size_t bins, bool density=false, bool return_edges=false);

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
