
#ifndef TENSOR_PYBIND11_H
#define TENSOR_PYBIND11_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// =================================================================================================
// type caster: cppmat::tensor4 <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<cppmat::tensor4<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(cppmat::tensor4<T>, _("cppmat::tensor4<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims != 4)
          return false;

        size_t nd = buf.shape()[0];

        for ( size_t i=0 ; i<buf.ndim() ; i++ )
          if ( buf.shape()[i] != nd )
            return false;

        value = cppmat::tensor4<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const cppmat::tensor4<T>& src,
        py::return_value_policy policy, py::handle parent)
      {
        std::vector<size_t> shape(4);
        for ( size_t i=0; i<4; ++i ) shape[i] = src.ndim();

        py::array a(std::move(shape),std::move(src.strides(true)),src.data());

        return a.release();
      }
  };
}} // namespace pybind11::detail

// =================================================================================================
// type caster: cppmat::tensor2 <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<cppmat::tensor2<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(cppmat::tensor2<T>, _("cppmat::tensor2<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims != 2)
          return false;

        size_t nd = buf.shape()[0];

        for ( size_t i=0 ; i<buf.ndim() ; i++ )
          if ( buf.shape()[i] != nd )
            return false;

        value = cppmat::tensor2<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const cppmat::tensor2<T>& src,
        py::return_value_policy policy, py::handle parent)
      {
        std::vector<size_t> shape(2);
        for ( size_t i=0; i<2; ++i ) shape[i] = src.ndim();

        py::array a(std::move(shape),std::move(src.strides(true)),src.data());

        return a.release();
      }
  };
}} // namespace pybind11::detail

// =================================================================================================
// type caster: cppmat::vector <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<cppmat::vector<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(cppmat::vector<T>, _("cppmat::vector<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims != 1)
          return false;

        size_t nd = buf.shape()[0];

        for ( size_t i=0 ; i<buf.ndim() ; i++ )
          if ( buf.shape()[i] != nd )
            return false;

        value = cppmat::vector<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const cppmat::vector<T>& src,
        py::return_value_policy policy, py::handle parent)
      {
        std::vector<size_t> shape(1);
        for ( size_t i=0; i<1; ++i ) shape[i] = src.ndim();

        py::array a(std::move(shape),std::move(src.strides(true)),src.data());

        return a.release();
      }
  };
}} // namespace pybind11::detail

#endif
