
.. _compile:

*********
Compiling
*********

Introduction
============

As described, this module is header only. So one just has to ``#include <cppmat/matrix.h>`` and/or ``#include <cppmat/tensor.h>`` somewhere in the source code, and to tell the compiler where the header-files are. For the latter, several ways are described below.

Before proceeding, a words about optimization. Of course one should use optimization when compiling the release of the code (``-O2`` or ``-O3``). But it is also a good idea to switch of the assertions in the code (mostly checks on size) that facilitate easy debugging, but do cost time. Therefore, include the flag ``-DNDEBUG``. Note that this is all C++ standard. I.e. it should be no surprise, and it always a good idea to do.

.. _compile_pkg-config:

pkg-config
==========

To simplify matters greatly one can use ``pkg-config`` to keep track of the location of the header files. To that matter on has to:

1. Copy the file ``cppmat.pc.in`` to ``cppmat.pc`` to some location that can be found by ``pkg_config`` (for example by adding ``export PKG_CONFIG_PATH=/path/to/cppmat.pc:$PKG_CONFIG_PATH`` to the ``.bashrc``).

2. Modify the line ``prefix=@CMAKE_INSTALL_PREFIX@`` to ``prefix=/path/to/cppmat``.

GNU / Clang
===========

Add the following compiler's arguments:

.. code-block:: bash

  -I${PATH_TO_CPPMAT}/include -std=c++11

(or ``-std=c++14``).

If :ref:`compile_pkg-config` is configured, one can also use

.. code-block:: bash

  `pkg-config --cflags cppmat`

instead of the first argument.


.. note:: **Tip**

  If you want to avoid separately including the header files using a compiler flag, ``git submodule`` is a nice way to go:

  1.  Include this module as a submodule using ``git submodule add https://github.com/tdegeus/cppmat.git``.

  2.  Replace the first line of this example by ``#include "cppmat/include/cppmat/matrix.h"``.

      *If you decide to manually copy the header file, you might need to modify this relative path to your liking.*

cmake
=====

Add the following to your ``CMakeLists.txt``:

.. code-block:: cmake

  find_package(PkgConfig)
  pkg_check_modules(CPPMAT REQUIRED cppmat)
  include_directories(${CPPMAT_INCLUDE_DIRS})

Obviously one should configure :ref:`compile_pkg-config` for this to work.
