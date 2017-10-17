/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/cppmat

================================================================================================= */

#ifndef CPPMAT_MATRIX_PYBIND11_H
#define CPPMAT_MATRIX_PYBIND11_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "matrix.h"
#include "macros.h"

namespace py = pybind11;

// =================================================================================================
// type caster: cppmat::matrix <-> NumPy-array
// =================================================================================================

namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<cppmat::matrix<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(cppmat::matrix<T>, _("cppmat::matrix<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert)
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims < 1)
          return false;

        std::vector<size_t> shape(buf.ndim());

        for ( ssize_t i=0 ; i<buf.ndim() ; i++ )
          shape[i] = buf.shape()[i];

        value = cppmat::matrix<T>(shape,buf.data());

        return true;
      }

      // Conversion part 2 (C++ -> Python)
      static py::handle cast(const cppmat::matrix<T>& src,
        py::return_value_policy policy, py::handle parent)
      {
        py::array a(std::move(src.shape()),std::move(src.strides(true)),src.data());

        return a.release();
      }
  };
}} // namespace pybind11::detail

#endif
