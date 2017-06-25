
#include <cppmat/tensor.h>
#include <cppmat/pybind11_tensor.h>

// =============================================================================
// double tensor contraction 4-d : 2-d -> 2-d                         (pure C++)
// =============================================================================

template <class T>
cppmat::tensor2<T> ddot42 ( cppmat::tensor4<T> &A, cppmat::tensor2<T> &B )
{
  return A.ddot(B);
}

template cppmat::tensor2<int>    ddot42<int>    (cppmat::tensor4<int>    &, cppmat::tensor2<int>    &);
template cppmat::tensor2<double> ddot42<double> (cppmat::tensor4<double> &, cppmat::tensor2<double> &);


// =============================================================================
// create Python module
// =============================================================================

PYBIND11_PLUGIN(tensorlib) {

py::module m("tensorlib", "Tensor library");

m.def("ddot42",py::overload_cast<cppmat::tensor4<int>   &,cppmat::tensor2<int>   &>(&ddot42<int>   ),"Double-dot product",py::arg("A"),py::arg("B"));
m.def("ddot42",py::overload_cast<cppmat::tensor4<double>&,cppmat::tensor2<double>&>(&ddot42<double>),"Double-dot product",py::arg("A"),py::arg("B"));

return m.ptr();

} // PYBIND11_PLUGIN
