
.. _compile:

*********
Compiling
*********

Introduction
============

As described, this module is header only. So one just has to ``#include <cppmat/matrix.h>`` and/or ``#include <cppmat/tensor.h>`` somewhere in the source code, and to tell the compiler where the header-files are. For the latter, several ways are described below.

Before proceeding, a words about optimization. Of course one should use optimization when compiling the release of the code (``-O2`` or ``-O3``). But it is also a good idea to switch of the assertions in the code (mostly checks on size) that facilitate easy debugging, but do cost time. Therefore, include the flag ``-DNDEBUG``. Note that this is all C++ standard. I.e. it should be no surprise, and it always a good idea to do.

GNU / Clang
===========

Add the following compiler's arguments:

.. code-block:: bash

  -I${PATH_TO_CPPMAT}/include -std=c++11

(or ``-std=c++14``).

.. note:: **Tip**

  If you want to avoid separately including the header files using a compiler flag, ``git submodule`` is a nice way to go:

  1.  Include this module as a submodule using ``git submodule add https://github.com/tdegeus/cppmat.git``.

  2.  Replace the first line of this example by ``#include "cppmat/include/cppmat/matrix.h"``.

      *If you decide to manually copy the header file, you might need to modify this relative path to your liking.*

  Or see :ref:`compile_automatic`. You can also combine the ``git submodule`` with any of the below compiling strategies.

.. _compile_automatic:

Automating build
================

Install
-------

To enable automatic build one should 'install' ``cppmat`` somewhere.

.. note::

  If you do not wish to use ``CMake``, or you want to do something custom. You can of course. Follow these steps:

  1.  Copy the file ``cppmat.pc.in`` to ``cppmat.pc`` to some location that can be found by ``pkg_config`` (for example by adding ``export PKG_CONFIG_PATH=/path/to/cppmat.pc:$PKG_CONFIG_PATH`` to the ``.bashrc``).

  2.  Modify the line ``prefix=@CMAKE_INSTALL_PREFIX@`` to ``prefix=/path/to/cppmat``.

  3.  Modify the line ``Version: @CPPMAT_VERSION_NUMBER@`` to reflect the correct release version.

Install system-wide (root)
^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Make a temporary build directory. For example

    .. code-block:: bash

      $ cd /path/to/cppmat
      $ mkdir build
      $ cd build

2.  'Build' ``cppmat``

    .. code-block:: bash

      $ cmake ..
      $ make install

Install in custom location (user)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Make a temporary build directory. For example

    .. code-block:: bash

      $ cd /path/to/cppmat
      $ mkdir build
      $ cd build

2.  'Build' ``cppmat``, to install it in a custom location

    .. code-block:: bash

      $ mkdir /custom/install/path
      $ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/custom/install/path
      $ make install

3.  Add the following path to your ``~/.bashrc`` (or ``~/.zshrc``):

    .. code-block:: bash

      export PKG_CONFIG_PATH=/custom/install/path/share/pkgconfig:$PKG_CONFIG_PATH

pkg-config
----------

Instead of ``-I...`` one can now use

.. code-block:: bash

  `pkg-config --cflags cppmat` -std=c++11

to compile in a single command.

cmake
-----

Add the following to your ``CMakeLists.txt``:

.. code-block:: cmake

  find_package(PkgConfig)
  pkg_check_modules(CPPMAT REQUIRED cppmat)
  include_directories(${CPPMAT_INCLUDE_DIRS})

