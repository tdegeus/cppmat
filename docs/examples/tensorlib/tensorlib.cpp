
#include "cpptensor.h"
#include "py_cpptensor.h"

// =============================================================================
// double tensor contraction 4-d : 2-d -> 2-d                         (pure C++)
// =============================================================================

template <class T>
tensor::tensor2<T> ddot42 ( tensor::tensor4<T> &A, tensor::tensor2<T> &B )
{
  return A.ddot(B);
}

template tensor::tensor2<int>    ddot42<int>    (tensor::tensor4<int>    &, tensor::tensor2<int>    &);
template tensor::tensor2<double> ddot42<double> (tensor::tensor4<double> &, tensor::tensor2<double> &);


// =============================================================================
// create Python module
// =============================================================================

PYBIND11_PLUGIN(tensorlib) {

py::module m("tensorlib", "Tensor library");

m.def("ddot42",py::overload_cast<tensor::tensor4<int>   &,tensor::tensor2<int>   &>(&ddot42<int>   ),"Double-dot product",py::arg("A"),py::arg("B"));
m.def("ddot42",py::overload_cast<tensor::tensor4<double>&,tensor::tensor2<double>&>(&ddot42<double>),"Double-dot product",py::arg("A"),py::arg("B"));

return m.ptr();

} // PYBIND11_PLUGIN
