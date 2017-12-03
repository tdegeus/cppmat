/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_TENSOR_PYBIND11_H
#define CPPMAT_TENSOR_PYBIND11_H

#include "pybind11.h"

namespace py = pybind11;

namespace pybind11 {
namespace detail {

// =================================================================================================
// type caster: cppmat::cartesian::tensor4 <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian::tensor4<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::tensor4<T>, _("cppmat::cartesian::tensor4<T>"));

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

    // - read number of dimensions (shape in each direction)
    ssize_t nd = buf.shape()[0];

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::tensor4<T>(nd);

    // - copy data
    std::copy(buf.data(), buf.data()+nd*nd*nd*nd, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::tensor4<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian::tensor2 <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian::tensor2<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::tensor2<T>, _("cppmat::cartesian::tensor2<T>"));

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

    // - read number of dimensions (shape in each direction)
    ssize_t nd = buf.shape()[0];

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::tensor2<T>(nd);

    // - copy data
    std::copy(buf.data(), buf.data()+nd*nd, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::tensor2<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - create Python variable (all variables are copied)
    py::array a(std::move(src.shape()), std::move(src.strides(true)), src.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian::tensor2s <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian::tensor2s<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::tensor2s<T>, _("cppmat::cartesian::tensor2s<T>"));

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

    // - read number of dimensions (shape in each direction)
    ssize_t nd = buf.shape()[0];

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::tensor2s<T>(nd);

    // - check for symmetry
    #ifndef NDEBUG
    for ( size_t i = 0 ; i < nd ; ++i )
      for ( size_t j = i+1 ; j < nd ; ++j )
        assert( buf.data()[i*nd+j] == buf.data()[j*nd+i] );
    #endif

    // - copy from input (ignores lower diagonal terms)
    for ( size_t i = 0 ; i < nd ; ++i )
      for ( size_t j = i ; j < nd ; ++j )
        value[i*nd-(i-1)*i/2+j-i] = buf.data()[i*nd+j];

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::tensor2s<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - convert to dense tensor
    cppmat::cartesian::tensor2<T> tmp = src;

    // - create Python variable (all variables are copied)
    py::array a(std::move(tmp.shape()), std::move(tmp.strides(true)), tmp.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian::tensor2d <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian::tensor2d<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::tensor2d<T>, _("cppmat::cartesian::tensor2d<T>"));

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

    // - read number of dimensions (shape in each direction)
    ssize_t nd = buf.shape()[0];

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::tensor2d<T>(nd);

    // - check the input to be diagonal
    #ifndef NDEBUG
    for ( size_t i = 0 ; i < nd ; ++i )
      for ( size_t j = 0 ; j < nd ; ++j )
        if ( i !=j )
          assert( !buf.data()[i*nd+j] );
    #endif

    // - copy from input (ignores off-diagonal terms)
    for ( size_t i = 0 ; i < nd ; ++i )
      value[i] = buf.data()[i*nd+i];

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::tensor2d<T>& src, py::return_value_policy policy, py::handle parent
  )
  {
    // - convert to dense tensor
    cppmat::cartesian::tensor2<T> tmp = src;

    // - create Python variable (all variables are copied)
    py::array a(std::move(tmp.shape()), std::move(tmp.strides(true)), tmp.data());

    // - release variable to Python
    return a.release();
  }
};

// =================================================================================================
// type caster: cppmat::cartesian::vector <-> NumPy-array
// =================================================================================================

template <typename T> struct type_caster<cppmat::cartesian::vector<T>>
{
public:

  PYBIND11_TYPE_CASTER(cppmat::cartesian::vector<T>, _("cppmat::cartesian::vector<T>"));

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

    // - read number of dimensions (shape in each direction)
    ssize_t nd = buf.shape()[0];

    // - the shape in each direction should be equal ( == nd )
    for ( ssize_t i = 0 ; i < rank ; ++i )
      if ( buf.shape()[i] != nd )
        return false;

    // - all checks passed : create the proper C++ variable
    value = cppmat::cartesian::vector<T>(nd);

    // - copy data
    std::copy(buf.data(), buf.data()+nd, value.data());

    // - signal successful variable creation
    return true;
  }

  // C++ -> Python
  // -------------

  static py::handle cast(
    const cppmat::cartesian::vector<T>& src, py::return_value_policy policy, py::handle parent
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
