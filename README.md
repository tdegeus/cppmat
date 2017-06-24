# cppmat

Header-only module that provides a C++ matrix class for n-d matrices.

```cpp
#include "cppmat.h"

int main()
{
    mat::matrix<double> A({10,10,10});
    // ...
    return 0;
}
```

Really, that's it! Not special compiler statements (except `-std=c++11` or `-std=c++14`), no libraries, nothing. 

# ccptensor

Header-only module that provides C++ classes for 4th- and 2nd order tensors and vectors (which are essentially vectors, but with special methods).

```cpp
#include "cpptensor"

int main()
{
    tensor::tensor4<double> tensor::identity4(3);
    // ...
    return 0;
}
```

# Python interface

This library includes a header that provides an interface to [pybind11](https://github.com/pybind/pybind11) such that an interface to NumPy arrays is automatically provided when including a function with a `mat::matrix`. One has to only `#include py_cppmat.h`. Likewise for functions with `tensor::tensor4`, `tensor::tensor2`, or `tensor::vector` one includes  `#include py_cpptensor.h`.

An example is provided in `docs/examples/tensorlib`. This example includes two forms of building:

1.  `CMakeList.txt` for building using `cmake` (`cmake .` and then `make`). For this to work, `pybind11` must be 'installed' on the system. Alternatively one include `pybind11` as a subfolder (for example using `git submodule add https://github.com/pybind/pybind11.git`). In that case, replace `find_package(pybind11 REQUIRED)` by `add_subdirectory(pybind11)` in `CMakeList.txt`.

2.  `setup.py` for building using `python` (`python3 setup.py build` and then `python3 setup.py install`). Using this option `python` will take care of the `pybind11` dependency.




