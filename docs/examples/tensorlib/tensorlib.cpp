
#include "cppmat.h"
#include "py_cppmat.h"

// =============================================================================
// double tensor contraction 4-d : 2-d -> 2-d                         (pure C++)
// =============================================================================

template <class T>
mat::matrix<T> ddot42 ( mat::matrix<T> &A, mat::matrix<T> &B )
{
  if ( A.ndim()!=4 || B.ndim()!=2 )
    throw std::length_error("A should be 4D, B should be 2D");
  if ( A.shape()[2]!=B.shape()[0] || A.shape()[3]!=B.shape()[1] )
    throw std::length_error("A and B inconsistent");

  mat::matrix<T> C({A.shape()[0],B.shape()[1]});
  C.zeros();

  for ( size_t i=0 ; i<A.shape()[0] ; ++i )
    for ( size_t j=0 ; j<A.shape()[1] ; ++j )
      for ( size_t k=0 ; k<A.shape()[2] ; ++k )
        for ( size_t l=0 ; l<A.shape()[3] ; ++l )
          C(i,j) += A(i,j,k,l)*B(l,k);

  return C;
}

template mat::matrix<int>    ddot42<int>    (mat::matrix<int>    &, mat::matrix<int>    &);
template mat::matrix<double> ddot42<double> (mat::matrix<double> &, mat::matrix<double> &);


// =============================================================================
// create Python module
// =============================================================================

PYBIND11_PLUGIN(tensorlib) {

py::module m("tensorlib", "Tensor library");

m.def("ddot42",py::overload_cast<mat::matrix<int>   &,mat::matrix<int>   &>(&ddot42<int>   ),"Double-dot product",py::arg("A"),py::arg("B"));
m.def("ddot42",py::overload_cast<mat::matrix<double>&,mat::matrix<double>&>(&ddot42<double>),"Double-dot product",py::arg("A"),py::arg("B"));

return m.ptr();

} // PYBIND11_PLUGIN
