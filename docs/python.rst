
.. _python:

****************
Python interface
****************

This library includes provides an interface to `pybind11 <https://github.com/pybind/pybind11>`_ such that an interface to NumPy arrays is automatically provided when including a function with any of the cppmat classes:

+---------------------------------------+---------------------------+
| **cppmat class**                      | **Rank of NumPy-array**   |
+=======================================+===========================+
| ``cppmat::array``                     | n                         |
+---------------------------------------+---------------------------+
| ``cppmat::matrix``                    | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::vector``                    | 1                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::array``               | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::matrix``              | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::vector``              | 1                         |
+---------------------------------------+---------------------------+
| ``cppmat::cartesian::tensor4``        | 4                         |
+---------------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2``        | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2s``       | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2d``       | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::cartesian::vector``         | 1                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::cartesian::tensor4``  | 4                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::cartesian::tensor2``  | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::cartesian::tensor2s`` | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::cartesian::tensor2d`` | 2                         |
+---------------------------------------+---------------------------+
| ``cppmat::tiny::cartesian::vector``   | 1                         |
+---------------------------------------+---------------------------+

.. warning::

  On the Python side all the matrices (``cppmat::matrix``, ``cppmat::symmetric::matrix``, and ``cppmat::diagonal::matrix``) and 2nd-order tensors (``cppmat::cartesian::tensor2``, ``cppmat::cartesian::tensor2s``, and ``cppmat::cartesian::tensor2d``) are all square matrices (rank 2 NumPy arrays). This means that when a function that has ``cppmat::symmetric::matrix`` or ``cppmat::cartesian::tensor2s`` as argument, the upper-diagonal part is read; while when it has an argument ``cppmat::diagonal::matrix`` or ``cppmat::cartesian::tensor2d`` only the diagonal is considered.

  **This requires extra attention as information might be lost. To optimize for speed and flexibility no checks are performed in the release libraries derived from cppmat!**

  You can ask cppmat to check for this, by omitting the ``-DNDEBUG`` compiler flag (this enables several assertions, so it may cost you some efficiency).

  (The same holds for the classes under ``cppmat::tiny::``.)

To use this feature one has to:

.. code-block:: cpp

  #include <cppmat/pybind11.h>

Building
========

[:download:`tensorlib.zip <./examples/tensorlib.zip>`]

Building is demonstrated based on the 'tensorlib' example.

CMake
-----

[:download:`CMakeLists.txt <./examples/tensorlib/CMakeLists.txt>`]

.. literalinclude:: examples/tensorlib/CMakeLists.txt
   :language: cmake

Build using

.. code-block:: bash

  cd /path/to/tempdir
  cmake /path/to/tensorlib
  make

For this to work, *pybind11* must be 'installed' on the system.

.. tip::

  Alternatively you can include *pybind11* as a sub-folder (for example using ``git submodule add https://github.com/pybind/pybind11.git``). In that case, replace ``find_package(pybind11 REQUIRED)`` by ``add_subdirectory(pybind11)`` in ``CMakeLists.txt``.

.. tip::

  To link to external libraries, include at the end of your ``CMakeLists.txt``

  .. code-block:: cmake

    target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_LIBS})

setup.py
--------

[:download:`setup.py <./examples/tensorlib/setup.py>`]

.. literalinclude:: examples/tensorlib/setup.py
   :language: python

As shown in this example building using the ``setup.py`` can be simplified using some routines in the ``cppmat`` python module. These routines have been taken from `pybind <https://github.com/pybind>`_, most notably from Sylvain Corlay and Dean Moldovan. They are merely attached to this module to simplify your life.

Build using

.. code-block:: bash

  python3 setup.py build
  python3 setup.py install

.. tip::

  Replace the executable with your favorite Python version, e.g. with ``python``.

CMake & setup.py
----------------

CMake can be called from the ``setup.py`` to take advantage of both. In that case the ``setup.py`` would be simply

.. code-block:: python

  import setuptools, cppmat

  setuptools.setup(
    name             = 'tensorlib',
    version          = '0.0.1',
    author           = 'Tom de Geus',
    author_email     = 'email@address.com',
    description      = 'Description',
    long_description = '',
    ext_modules      = [cppmat.CMakeExtension('tensorlib')],
    cmdclass         = dict(build_ext=cppmat.CMakeBuild),
    zip_safe         = False,
  )



