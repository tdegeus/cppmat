
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "cpptensor.h"

namespace py = pybind11;

// =================================================================================================
// type caster: tensor::tensor4 <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<tensor::tensor4<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(tensor::tensor4<T>, _("tensor::tensor4<T>"));

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

        value = tensor::tensor4<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const tensor::tensor4<T>& src,
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
// type caster: tensor::tensor2 <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<tensor::tensor2<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(tensor::tensor2<T>, _("tensor::tensor2<T>"));

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

        value = tensor::tensor2<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const tensor::tensor2<T>& src,
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
// type caster: tensor::vector <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<tensor::vector<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(tensor::vector<T>, _("tensor::vector<T>"));

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

        value = tensor::vector<T>(nd,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const tensor::vector<T>& src,
        py::return_value_policy policy, py::handle parent)
      {
        std::vector<size_t> shape(1);
        for ( size_t i=0; i<1; ++i ) shape[i] = src.ndim();

        py::array a(std::move(shape),std::move(src.strides(true)),src.data());

        return a.release();
      }
  };
}} // namespace pybind11::detail
