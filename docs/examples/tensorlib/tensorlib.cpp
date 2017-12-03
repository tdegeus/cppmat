
#include <cppmat/cppmat.h>
#include <cppmat/pybind11.h>

// =================================================================================================
// double tensor contraction 4-d : 2-d -> 2-d                                             (pure C++)
// =================================================================================================

namespace cm = cppmat::cartesian;

template <class T>
cm::tensor2<T> ddot42(cm::tensor4<T> &A, cm::tensor2<T> &B)
{
  return A.ddot(B);
}

template cm::tensor2<int>    ddot42<int>   (cm::tensor4<int>   &, cm::tensor2<int>   &);
template cm::tensor2<double> ddot42<double>(cm::tensor4<double>&, cm::tensor2<double>&);

// -------------------------------------------------------------------------------------------------

cm::tensor2s<double> symmetric(cm::tensor2<double> &B)
{
  return B.cast<cm::tensor2s<double>>();
}

// -------------------------------------------------------------------------------------------------

cm::tensor2s<double> diagonal(cm::tensor2<double> &B)
{
  return B.cast<cm::tensor2d<double>>();
}

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

m.def("symmetric",&symmetric,
  "Symmetric part of the tensor",
  py::arg("A")
);

m.def("diagonal",&diagonal,
  "Diagonal part of the tensor",
  py::arg("A")
);

} // PYBIND11_MODULE
