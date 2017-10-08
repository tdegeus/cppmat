
#include <cppmat/tensor.h>
#include <cppmat/pybind11_tensor.h>

// =================================================================================================
// double tensor contraction 4-d : 2-d -> 2-d                                             (pure C++)
// =================================================================================================

namespace cm = cppmat::cartesian;

template <class T>
cm::tensor2<T> ddot42 ( cm::tensor4<T> &A, cm::tensor2<T> &B )
{
  return A.ddot(B);
}

template cm::tensor2<int>    ddot42<int>   (cm::tensor4<int>   &, cm::tensor2<int>   &);
template cm::tensor2<double> ddot42<double>(cm::tensor4<double>&, cm::tensor2<double>&);

// -------------------------------------------------------------------------------------------------

template <class T>
cm::tensor2s<T> symmetric ( cm::tensor2<T> &B )
{
  cm::tensor2s<T> C = B.astensor2s();
  return C;
}

template cm::tensor2s<int>    symmetric<int>    (cm::tensor2<int>    &);
template cm::tensor2s<double> symmetric<double> (cm::tensor2<double> &);

// -------------------------------------------------------------------------------------------------

template <class T>
cm::tensor2d<T> diagonal ( cm::tensor2<T> &B )
{
  cm::tensor2d<T> C = B.astensor2d();
  return C;
}

template cm::tensor2d<int>    diagonal<int>    (cm::tensor2<int>    &);
template cm::tensor2d<double> diagonal<double> (cm::tensor2<double> &);

// =================================================================================================
// create Python module                                                                   (pybind11)
// =================================================================================================

PYBIND11_MODULE(tensorlib, m)
{

m.doc() = "Tensor library";

m.def("ddot42",py::overload_cast<cm::tensor4<int>&,cm::tensor2<int>&>(&ddot42<int>),
  "Double-dot product",
  py::arg("A"),
  py::arg("B")
);

m.def("ddot42",py::overload_cast<cm::tensor4<double>&,cm::tensor2<double>&>(&ddot42<double>),
  "Double-dot product",
  py::arg("A"),
  py::arg("B")
);


m.def("symmetric",py::overload_cast<cm::tensor2<int>&>(&symmetric<int>),
  "Symmetric part of the tensor",
  py::arg("A")
);

m.def("symmetric",py::overload_cast<cm::tensor2<double>&>(&symmetric<double>),
  "Symmetric part of the tensor",
  py::arg("A")
);


m.def("diagonal",py::overload_cast<cm::tensor2<int>&>(&diagonal<int>),
  "Diagonal part of the tensor",
  py::arg("A")
);

m.def("diagonal",py::overload_cast<cm::tensor2<double>&>(&diagonal<double>),
  "Diagonal part of the tensor",
  py::arg("A")
);

} // PYBIND11_MODULE
