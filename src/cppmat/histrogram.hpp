/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_HISTOGRAM_HPP
#define CPPMAT_HISTOGRAM_HPP

// -------------------------------------------------------------------------------------------------

#include "cppmat.h"

// -------------------------------------------------------------------------------------------------

namespace cppmat {

// =================================================================================================

template<typename X>
std::tuple<std::vector<double>, std::vector<double>> histogram(
  const std::vector<X> &data, size_t bins, bool density, bool return_edges)
{
  // alias
  double Bins = static_cast<double>(bins);
  double N    = static_cast<double>(data.size());

  // get domain, and discretization size
  double min = *std::min_element(data.begin(),data.end());
  double max = *std::max_element(data.begin(),data.end());
  double h   = (max - min) / Bins;

  // histogram
  // - zero-initialize
  std::vector<size_t> count(bins, 0);
  // - fill
  for ( auto &i : data )
    count[std::min(static_cast<size_t>((static_cast<double>(i)-min)/h),bins-1)]++;

  // type-cast
  std::vector<double> countD(count.begin(), count.end());

  // convert to density: set the integral to one
  if ( density )
    for ( auto &i : countD )
      i /= ( h * N );

  // return output
  if ( return_edges ) return std::make_tuple(countD, cppmat::linspace(min     , max     , bins+1));
  else                return std::make_tuple(countD, cppmat::linspace(min+h/2., max-h/2., bins  ));
}

// -------------------------------------------------------------------------------------------------

template<typename X>
std::tuple<std::vector<double>, std::vector<double>> histogram_uniform(
  const std::vector<X> &data, size_t bins, bool density, bool return_edges)
{
  // alias
  double Bins = static_cast<double>(bins);
  double N    = static_cast<double>(data.size());

  // histogram
  // - all bins equal
  std::vector<size_t> count(bins, std::floor(N/Bins));
  // - remainder
  size_t diff = data.size() - bins * static_cast<size_t>(count[0]);
  // - locations to distribute remainder
  std::vector<size_t> idx = cppmat::linspace<size_t>(0, bins-1, diff);
  // - distribute remainder
  for ( auto &i : idx ) count[i]++;

  // sorted data
  // - allocate
  std::vector<X> sorted = data;
  // - sort
  std::sort(sorted.begin(), sorted.end());

  // edges
  // - allocate
  std::vector<X> edges(bins+1);
  // - counter
  size_t j = 0;
  // - compute
  for ( size_t i = 0 ; i < bins ; ++i ) {
    edges[i] = sorted[j];
    j += count[i];
  }
  // - last edge
  edges[bins] = sorted[sorted.size()-1];

  // type-cast
  std::vector<double> countD(count.begin(), count.end());
  std::vector<double> edgesD(edges.begin(), edges.end());

  // convert to density: set the integral to one
  if ( density )
  {
    // - zero-initialize integral
    double I = 0.0;
    // - compute integral
    for ( size_t i = 0 ; i < bins ; ++i )
      I += ( edgesD[i+1] - edgesD[i] ) * countD[i];
    // - normalize
    for ( auto &i : countD )
      i /= I;
  }

  // return edges
  if ( return_edges ) return std::make_tuple(countD, edgesD);

  // return mid-points
  // - allocate
  std::vector<double> mid(bins);
  // - compute
  for ( size_t i = 0 ; i < bins ; ++i )
    mid[i] = ( edges[i+1] - edges[i] ) / 2.;
  // - return
  return std::make_tuple(countD, mid);
}

// =================================================================================================

} // namespace ...

// -------------------------------------------------------------------------------------------------

#endif
