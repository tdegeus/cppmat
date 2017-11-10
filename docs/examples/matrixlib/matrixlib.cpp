
#include <cppmat/cppmat.h>
#include <cppmat/pybind11.h>

// =================================================================================================
// sqrt of every entry                                                                    (pure C++)
// =================================================================================================

cppmat::matrix<double> squareRoot ( cppmat::matrix<double> &A )
{
  cppmat::matrix<double> C = A;

  for ( auto &i : C )
    i = std::sqrt(i);

  return C;
}

// =================================================================================================
// create Python module                                                                   (pybind11)
// =================================================================================================

PYBIND11_MODULE(matrixlib, m)
{

m.doc() = "Matrix library";

m.def("squareRoot",&squareRoot,
  "Square root of every entry",
  py::arg("A")
);

} // PYBIND11_MODULE
