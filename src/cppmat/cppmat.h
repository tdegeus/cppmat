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
#define CPPMAT_MAJOR_VERSION 4
#define CPPMAT_MINOR_VERSION 6

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

#include "matrix.h"
#include "matrix2.h"
#include "vector.h"
#include "periodic_matrix.h"
#include "periodic_matrix2.h"
#include "periodic_vector.h"
#include "tiny_matrix2.h"
#include "tiny_vector.h"
#include "view_matrix2.h"
#include "view_vector.h"
#include "tensor.h"
#include "tensor2.h"
#include "tensor3.h"
#include "view_tensor2.h"
#include "view_tensor3.h"
#include "stl.h"

#include "matrix.cpp"
#include "matrix2.cpp"
#include "vector.cpp"
#include "periodic_matrix.cpp"
#include "periodic_matrix2.cpp"
#include "periodic_vector.cpp"
#include "tiny_matrix2.cpp"
#include "tiny_vector.cpp"
#include "view_matrix2.cpp"
#include "view_vector.cpp"
#include "tensor.cpp"
#include "tensor2.cpp"
#include "tensor3.cpp"
#include "view_tensor2.cpp"
#include "view_tensor3.cpp"
#include "stl.cpp"

// =================================================================================================

#endif

