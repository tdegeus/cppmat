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
#include <random>
#include <ctime>
#include <iso646.h> // to fix a Microsoft Visual Studio error on "and" and "or"

// =================================================================================================

#define CPPMAT_WORLD_VERSION 1
#define CPPMAT_MAJOR_VERSION 0
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

  template<class X> class array;
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
namespace cartesian {

  template<class X> class tensor4;
  template<class X> class tensor2;
  template<class X> class tensor2s;
  template<class X> class tensor2d;
  template<class X> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {

  template<class X, size_t RANK, size_t I, size_t J=1, size_t K=1, size_t L=1, size_t M=1, size_t N=1> class array;
  template<class X, size_t M, size_t N> class matrix;
  template<class X, size_t M> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace symmetric {

  template<class X, size_t M, size_t N=M> class matrix;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace diagonal {

  template<class X, size_t M, size_t N=M> class matrix;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace tiny {
namespace cartesian {

  template<class X, size_t ND> class tensor4;
  template<class X, size_t ND> class tensor2;
  template<class X, size_t ND> class tensor2s;
  template<class X, size_t ND> class tensor2d;
  template<class X, size_t ND> class vector;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {

  template<class X, size_t RANK, size_t I, size_t J=1, size_t K=1, size_t L=1, size_t M=1, size_t N=1> class array;
  template<class X, size_t M, size_t N> class matrix;
  template<class X, size_t M> class vector;

}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace symmetric {

  template<class X, size_t M, size_t N=M> class matrix;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace diagonal {

  template<class X, size_t M, size_t N=M> class matrix;

}}}

// -------------------------------------------------------------------------------------------------

namespace cppmat {
namespace view {
namespace cartesian {

  template<class X, size_t ND> class tensor4;
  template<class X, size_t ND> class tensor2;
  template<class X, size_t ND> class tensor2s;
  template<class X, size_t ND> class tensor2d;
  template<class X, size_t ND> class vector;

}}}

// =================================================================================================

#include "stl.h"
#include "private.h"

#include "var_regular_array.h"
#include "var_regular_matrix.h"
#include "var_regular_vector.h"
#include "var_symmetric_matrix.h"
#include "var_diagonal_matrix.h"
#include "var_misc_matrix.h"
#include "var_cartesian.h"
#include "var_cartesian_tensor4.h"
#include "var_cartesian_tensor2.h"
#include "var_cartesian_tensor2s.h"
#include "var_cartesian_tensor2d.h"
#include "var_cartesian_vector.h"

#include "fix_regular_array.h"
#include "fix_regular_matrix.h"
#include "fix_regular_vector.h"
#include "fix_symmetric_matrix.h"
#include "fix_diagonal_matrix.h"
#include "fix_misc_matrix.h"
#include "fix_cartesian.h"
#include "fix_cartesian_tensor4.h"
#include "fix_cartesian_tensor2.h"
#include "fix_cartesian_tensor2s.h"
#include "fix_cartesian_tensor2d.h"
#include "fix_cartesian_vector.h"

#include "map_regular_array.h"
#include "map_regular_matrix.h"
#include "map_regular_vector.h"
#include "map_symmetric_matrix.h"
#include "map_diagonal_matrix.h"
#include "map_cartesian_tensor4.h"
#include "map_cartesian_tensor2.h"
#include "map_cartesian_tensor2s.h"
#include "map_cartesian_tensor2d.h"
#include "map_cartesian_vector.h"

#include "stl.hpp"
#include "private.hpp"

#include "var_regular_array.hpp"
#include "var_regular_matrix.hpp"
#include "var_regular_vector.hpp"
#include "var_symmetric_matrix.hpp"
#include "var_diagonal_matrix.hpp"
#include "var_misc_matrix.hpp"
#include "var_cartesian.hpp"
#include "var_cartesian_tensor4.hpp"
#include "var_cartesian_tensor2.hpp"
#include "var_cartesian_tensor2s.hpp"
#include "var_cartesian_tensor2d.hpp"
#include "var_cartesian_vector.hpp"

#include "fix_regular_array.hpp"
#include "fix_regular_matrix.hpp"
#include "fix_regular_vector.hpp"
#include "fix_symmetric_matrix.hpp"
#include "fix_diagonal_matrix.hpp"
#include "fix_misc_matrix.hpp"
#include "fix_cartesian.hpp"
#include "fix_cartesian_2.hpp"
#include "fix_cartesian_3.hpp"
#include "fix_cartesian_tensor4.hpp"
#include "fix_cartesian_tensor2.hpp"
#include "fix_cartesian_tensor2s.hpp"
#include "fix_cartesian_tensor2d.hpp"
#include "fix_cartesian_vector.hpp"

#include "map_regular_array.hpp"
#include "map_regular_matrix.hpp"
#include "map_regular_vector.hpp"
#include "map_symmetric_matrix.hpp"
#include "map_diagonal_matrix.hpp"
#include "map_cartesian_tensor4.hpp"
#include "map_cartesian_tensor2.hpp"
#include "map_cartesian_tensor2s.hpp"
#include "map_cartesian_tensor2d.hpp"
#include "map_cartesian_vector.hpp"

// =================================================================================================

#endif

