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
#include <numeric>
#include <iso646.h> // to fix a Microsoft Visual Studio error on "and" and "or"

// =================================================================================================

#define CPPMAT_WORLD_VERSION 0
#define CPPMAT_MAJOR_VERSION 7
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

  template<class X> class Expandable;

}

// -------------------------------------------------------------------------------------------------

namespace cppmat {

  template<class X> class array;
  template<class X> class matrix4;
  template<class X> class matrix;
  template<class X> class vector;

}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace symmetric {

  template<class X> class matrix;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace diagonal {

  template<class X> class matrix;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace periodic {

  template<class X> class array;
  template<class X> class matrix;
  template<class X> class vector;

}}

// -------------------------------------------------------------------------------------------------

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
namespace cartesian2d {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}}

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

#include "stl.h"
#include "private.h"
#include "regular_array.h"
#include "regular_matrix.h"
#include "regular_vector.h"
#include "symmetric_matrix.h"
#include "diagonal_matrix.h"
#include "misc_matrix.h"

#include "stl.cpp"
#include "private.cpp"
#include "regular_array.cpp"
#include "regular_matrix.cpp"
#include "regular_vector.cpp"
#include "symmetric_matrix.cpp"
#include "diagonal_matrix.cpp"
#include "misc_matrix.cpp"

// =================================================================================================

#endif

