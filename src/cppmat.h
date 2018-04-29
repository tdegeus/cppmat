/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_H
#define CPPMAT_H

// =================================================================================================

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <iso646.h> // to fix a Microsoft Visual Studio error on "and" and "or"

// =================================================================================================

#define CPPMAT_WORLD_VERSION 0
#define CPPMAT_MAJOR_VERSION 6
#define CPPMAT_MINOR_VERSION 0

#define CPPMAT_VERSION_AT_LEAST(x,y,z) \
  (CPPMAT_WORLD_VERSION>x || (CPPMAT_WORLD_VERSION>=x && \
  (CPPMAT_MAJOR_VERSION>y || (CPPMAT_MAJOR_VERSION>=y && \
                              CPPMAT_MINOR_VERSION>=z))))

#define CPPMAT_VERSION(x,y,z) \
  (CPPMAT_WORLD_VERSION==x && \
   CPPMAT_MAJOR_VERSION==y && \
   CPPMAT_MINOR_VERSION==z)

// =================================================================================================

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// ====================================== forward declaration ======================================

namespace cppmat {
namespace cartesian {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian2d {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian2d {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace cartesian3d {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian3d {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}}

// =================================================================================================

#include "cppmat/matrix.h"
#include "cppmat/matrix2.h"
#include "cppmat/vector.h"
#include "cppmat/periodic_matrix.h"
#include "cppmat/periodic_matrix2.h"
#include "cppmat/periodic_vector.h"
#include "cppmat/tiny_matrix2.h"
#include "cppmat/tiny_vector.h"
#include "cppmat/view_matrix2.h"
#include "cppmat/view_vector.h"
#include "cppmat/tensor.h"
#include "cppmat/tensor2.h"
#include "cppmat/tensor3.h"
#include "cppmat/view_tensor2.h"
#include "cppmat/view_tensor3.h"
#include "cppmat/stl.h"

#include "cppmat/matrix.cpp"
#include "cppmat/matrix2.cpp"
#include "cppmat/vector.cpp"
#include "cppmat/periodic_matrix.cpp"
#include "cppmat/periodic_matrix2.cpp"
#include "cppmat/periodic_vector.cpp"
#include "cppmat/tiny_matrix2.cpp"
#include "cppmat/tiny_vector.cpp"
#include "cppmat/view_matrix2.cpp"
#include "cppmat/view_vector.cpp"
#include "cppmat/tensor.cpp"
#include "cppmat/tensor2.cpp"
#include "cppmat/tensor3.cpp"
#include "cppmat/view_tensor2.cpp"
#include "cppmat/view_tensor3.cpp"
#include "cppmat/stl.cpp"

// =================================================================================================

#endif

