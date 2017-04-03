# cppmat

Header-only module that provides a C++ matrix class for n-d matrices.

```cpp
#include "cppmat.h"

int main ( void )
{
    mat::matrix<double> A({10,10,10});
    // ...
    return 0;
}

```

Really, that's it! Not special compiler statements (except `-std=c++11` or `-std=c++14`), no libraries, nothing. 

# Python interface

This library includes a header that provides an interface to [pybind11](https://github.com/pybind/pybind11) such that an interface to NumPy arrays is automatically provided when including a function with a `mat::matrix`. One has to only `#include py_cppmat.h` (see *tensorlib example*).




