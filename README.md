
- [cppmat](#cppmat)
- [ccptensor](#ccptensor)
- [Compiling](#compiling)
- [Python interface](#python-interface)
- [Develop](#develop)

# cppmat

Header-only module that provides a C++ class for n-d matrices.

```cpp
#include <cppmat/matrix.h>

int main()
{
    cppmat::matrix<double> A({10,10,10});
    // ...
    return 0;
}
```

Really, that's it! You only need to [tell your compiler](#compiling) where the header files are (and to use a modern C++ standard). 

> If you want to avoid using compiler flags, `git submodule` is a nice way to go:
> 
>  1.  Include this module as a submodule using `git submodule add https://github.com/tdegeus/cppmat.git`.
>  2.  Replace the first line of this example by `#include "cppmat/include/cppmat/matrix.h"`.
>  
> *If you decide to manually copy the header file, you might need to modify this relative path to your liking.*

# ccptensor

Header-only module that provides C++ classes for 4th- and 2nd order tensors and vectors (the latter essentially coincide with `std::vector`, but with special methods).

```cpp
#include <cppmat/tensor.h>

int main()
{
    cppmat::tensor4<double> cppmat::identity4(3);
    // ...
    return 0;
}
```

See compilation remarks [above](#cppmat), and details [below](#compiling).

# Compiling

## pkg-config

To simplify matters greatly one can use `pkg-config` to keep track of the location of the header files. To that matter on has to:

1. Copy the file `cppmat.pc.in` to `cppmat.pc` to some location that can be found by `pkg_config` (for example by adding `export PKG_CONFIG_PATH=/path/to/cppmat.pc:$PKG_CONFIG_PATH` to the `.bashrc`). 
2. Modify the line `prefix=@CMAKE_INSTALL_PREFIX@` to `prefix=/path/to/cppmat`.

## GNU / Clang

Add the following compiler's arguments:

```bash
-I${PATH_TO_CPPMAT}/include -std=c++11
```

(or `-std=c++14`).

If `pkg-config` is configured on can also use

```bash
`pkg-config --cflags cppmat`
```

## cmake

Add the following to your `CMakeLists.txt`:

```cmake
find_package(PkgConfig)
pkg_check_modules(CPPMAT REQUIRED cppmat)
include_directories(${CPPMAT_INCLUDE_DIRS})
```

# Python interface

This library includes provides an interface to [pybind11](https://github.com/pybind/pybind11) such that an interface to NumPy arrays is automatically provided when including a function with a `cppmat::matrix`, `tensor::tensor4`, `tensor::tensor2`, or `tensor::vector`. To use this feature one has to include (either or both):

```cpp
#include <cppmat/pybind11_matrix.h>
#include <cppmat/pybind11_tensor.h>
```

An example is provided in `docs/examples/tensorlib`. This example includes two forms of building:

1.  `CMakeList.txt` for building using `cmake` (`cmake .` and then `make`). For this to work, `pybind11` must be 'installed' on the system. Alternatively you can include `pybind11` as a sub-folder (for example using `git submodule add https://github.com/pybind/pybind11.git`). In that case, replace `find_package(pybind11 REQUIRED)` by `add_subdirectory(pybind11)` in `CMakeList.txt`.

2.  `setup.py` for building using `python` (`python3 setup.py build` and then `python3 setup.py install`). Using this option `python` will take care of the `pybind11` and `cppmat` dependencies.

# Develop

## Python

The Python package of this module `cppmat/__init__.py` is essentially used to allow distribution of the header files that constitute this library through PyPi. In addition a small Python package `cppmat` is provided that allows easy `setup.py` formulations of derived packages. These features can also be used when one is just interested in using pybind11 and one does not intend to use `cppmat` itself.

## Create a new release

1.  Update the version numbers as follows:

    -   Modify `__version__` in `setup.py`.
    -   Modify `version` in `cppmat.pc.in`

2.  Upload the changes to GitHub and create a new release there (with the correct version number).

3.  Upload the package to PyPi:

    ```bash
    $ python3 setup.py bdist_wheel --universal
    $ twine upload dist/*
    ```

