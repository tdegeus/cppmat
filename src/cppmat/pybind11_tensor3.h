/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR3_PYBIND11_H
#define CPPMAT_TENSOR3_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::cartesian3d::tensor4 <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian3d::tensor4<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian3d::tensor4<T>, _("cppmat::cartesian3d::tensor4<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 4
    auto rank = buf.ndim();
    // - check
    if ( rank != 4 ) return false;

    // - read number of dimensions (shape in each direction) : should be exactly 3
    ssize_t nd = buf.shape()[0];
    // - check
    if ( nd != 3 ) return false;

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian2d::tensor4<T>();

    // - copy data
    std::copy(buf.data(), buf.data()+81, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian3d::tensor4<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian3d::tensor2 <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian3d::tensor2<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian3d::tensor2<T>, _("cppmat::cartesian3d::tensor2<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 2
    auto rank = buf.ndim();
    // - check
    if ( rank != 2 ) return false;

    // - read number of dimensions (shape in each direction) : should be exactly 3
    ssize_t nd = buf.shape()[0];
    // - check
    if ( nd != 3 ) return false;

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian3d::tensor2<T>();

    // - copy data
    std::copy(buf.data(), buf.data()+9, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian3d::tensor2<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian3d::tensor2s <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian3d::tensor2s<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian3d::tensor2s<T>, _("cppmat::cartesian3d::tensor2s<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 2
    auto rank = buf.ndim();
    // - check
    if ( rank != 2 ) return false;

    // - read number of dimensions (shape in each direction) : should be exactly 3
    ssize_t nd = buf.shape()[0];
    // - check
    if ( nd != 3 ) return false;

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian3d::tensor2s<T>();

    // - check for symmetry
    assert( buf[1] == buf[3] );
    assert( buf[2] == buf[6] );
    assert( buf[5] == buf[7] );

    // - copy from input (ignores lower diagonal terms)
    value[0] = buf[0];
    value[1] = buf[1];
    value[2] = buf[2];
    value[3] = buf[4];
    value[4] = buf[5];
    value[5] = buf[8];

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian3d::tensor2s<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - convert to dense tensor
    cppmat::cartesian3d::tensor2<T> tmp = src;

    // - create Python variable (all variables are copied)
    py::array a(std::move(tmp.shape()), std::move(tmp.strides(true)), tmp.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian3d::tensor2d <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian3d::tensor2d<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian3d::tensor2d<T>, _("cppmat::cartesian3d::tensor2d<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 2
    auto rank = buf.ndim();
    // - check
    if ( rank != 2 ) return false;

    // - read number of dimensions (shape in each direction) : should be exactly 3
    ssize_t nd = buf.shape()[0];
    // - check
    if ( nd != 3 ) return false;

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian3d::tensor2d<T>();

    // - check the input to be diagonal
    assert( ! buf[1] );
    assert( ! buf[2] );
    assert( ! buf[3] );
    assert( ! buf[5] );
    assert( ! buf[6] );
    assert( ! buf[7] );

    // - copy from input (ignores off-diagonal terms)
    value[0] = buf[0];
    value[1] = buf[4];
    value[2] = buf[8];

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian3d::tensor2d<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - convert to dense tensor
    cppmat::cartesian3d::tensor2<T> tmp = src;

    // - create Python variable (all variables are copied)
    py::array a(std::move(tmp.shape()), std::move(tmp.strides(true)), tmp.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian3d::vector <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian3d::vector<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian3d::vector<T>, _("cppmat::cartesian3d::vector<T>"));

  // Python -> C++
  // -------------

  bool load(py::handle src, bool convert)
  {
    // - basic pybind11 check
    if ( !convert && !py::array_t<T>::check_(src) ) return false;

    // - storage requirements : contiguous and row-major storage from NumPy
    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
    // - check
    if ( !buf ) return false;

    // - rank of the input array (number of indices) : should be exactly 1
    auto rank = buf.ndim();
    // - check
    if ( rank != 1 ) return false;

    // - read number of dimensions (shape in each direction) : should be exactly 3
    ssize_t nd = buf.shape()[0];
    // - check
    if ( nd != 3 ) return false;

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian3d::vector<T>();

    // - copy data
    std::copy(buf.data(), buf.data()+3, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian3d::vector<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================

}} // namespace pybind11::detail

#endif
