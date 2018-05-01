
.. _compile:

*********
Compiling
*********

Introduction
============

This module is header only. So one just has to ``#include <cppmat/cppmat.h>`` somewhere in the source code, and to tell the compiler where the header files are. For the latter, several ways are described below.

Before proceeding, a word about optimization. Of course one should use optimization when compiling the release of the code (``-O2`` or ``-O3``). But it is also a good idea to switch off the assertions in the code (mostly checks on size) that facilitate easy debugging, but do cost time. Therefore, include the flag ``-DNDEBUG``. Note that this is all C++ standard. I.e. it should be no surprise, and it is always a good idea to do.

Manual compiler flags
=====================

GNU / Clang
-----------

Add the following compiler's arguments:

.. code-block:: bash

  -I${PATH_TO_CPPMAT}/src -std=c++14

.. note:: **(Not recommended)**

  If you want to avoid separately including the header files using a compiler flag, ``git submodule`` is a nice way to go:

  1.  Include this module as a submodule using ``git submodule add https://github.com/tdegeus/cppmat.git``.

  2.  Replace the first line of this example by ``#include "cppmat/src/cppmat/cppmat.h"``.

      *If you decide to manually copy the header file, you might need to modify this relative path to your liking.*

  Or see :ref:`compile_automatic`. You can also combine the ``git submodule`` with any of the below compiling strategies.

.. _compile_automatic:

(Semi-)Automatic compiler flags
===============================

Install
-------

To enable (semi-)automatic build, one should 'install' ``cppmat`` somewhere.

Install systemwide (root)
^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Proceed to a (temporary) build directory. For example

    .. code-block:: bash

      $ cd /path/to/temp/build

2.  'Build' ``cppmat``

    .. code-block:: bash

      $ cmake /path/to/cppmat
      $ make install

.. note::

  One usually does not need any compiler arguments after following this protocol.

Install in custom location (user)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1.  Proceed to a (temporary) build directory. For example

    .. code-block:: bash

      $ cd /path/to/temp/build

2.  'Build' ``cppmat``, to install it in a custom location

    .. code-block:: bash

      $ mkdir /custom/install/path
      $ cmake /path/to/cppmat -DCMAKE_INSTALL_PREFIX:PATH=/custom/install/path
      $ make install

3.  Add the following path to your ``~/.bashrc`` (or ``~/.zshrc``):

    .. code-block:: bash

      export PKG_CONFIG_PATH=/custom/install/path/share/pkgconfig:$PKG_CONFIG_PATH
      export CPLUS_INCLUDE_PATH=$HOME/custom/install/path/include:$CPLUS_INCLUDE_PATH

.. note::

  One usually does not need any compiler arguments after following this protocol.

.. note:: **(Not recommended)**

  If you do not wish to use ``CMake`` for the installation, or you want to do something custom. You can, of course. Follow these steps:

  1.  Copy the file ``src/cppmat.pc.in`` to ``cppmat.pc`` to some location that can be found by ``pkg_config`` (for example by adding ``export PKG_CONFIG_PATH=/path/to/cppmat.pc:$PKG_CONFIG_PATH`` to the ``.bashrc``).

  2.  Modify the line ``prefix=@CMAKE_INSTALL_PREFIX@`` to ``prefix=/path/to/cppmat``.

  3.  Modify the line ``Cflags: -I${prefix}/@CPPMAT_INCLUDE_DIR@`` to ``Cflags: -I${prefix}/src``.

  4.  Modify the line ``Version: @CPPMAT_VERSION_NUMBER@`` to reflect the correct release version.

Compiler arguments from 'pkg-config'
------------------------------------

Should the compiler for some reason not be able to find the headers, instead of ``-I...`` one can now use

.. code-block:: bash

  `pkg-config --cflags cppmat` -std=c++14

as compiler argument.

Compiler arguments from 'cmake'
-------------------------------

Add the following to your ``CMakeLists.txt``:

.. code-block:: cmake

  set(CMAKE_CXX_STANDARD 14)

  find_package(PkgConfig)

  pkg_check_modules(CPPMAT REQUIRED cppmat)
  include_directories(${CPPMAT_INCLUDE_DIRS})

.. note::

  Except the C++ standard it should usually not be necessary to load cppmat explicitly, as it is installed in a location when the compiler can find it.
