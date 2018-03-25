
.. _python:

****************
Python interface
****************

This library includes provides an interface to `pybind11 <https://github.com/pybind/pybind11>`_ such that an interface to NumPy arrays is automatically provided when including a function with any of the cppmat classes:

+-----------------------------------+---------------------------+
| **cppmat class**                  | **Rank of NumPy-array**   |
+===================================+===========================+
| ``cppmat::matrix``                | n                         |
+-----------------------------------+---------------------------+
| ``cppmat::matrix2``               | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::vector``                | 1                         |
+-----------------------------------+---------------------------+
| ``cppmat::tiny::matrix2``         | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::tiny::vector``          | 1                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian::tensor4``    | 4                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2``    | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2s``   | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian::tensor2d``   | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian::vector``     | 1                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian2d::tensor4``  | 4                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian2d::tensor2``  | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian2d::tensor2s`` | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian2d::tensor2d`` | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian2d::vector``   | 1                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian3d::tensor4``  | 4                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian3d::tensor2``  | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian3d::tensor2s`` | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian3d::tensor2d`` | 2                         |
+-----------------------------------+---------------------------+
| ``cppmat::cartesian3d::vector``   | 1                         |
+-----------------------------------+---------------------------+


.. note:: **Warning**

  On the Python side all 2nd-order tensors (``cppmat::cartesian::tensor2``, ``cppmat::cartesian::tensor2s``, and ``cppmat::cartesian::tensor2d``) are all square matrices (rank 2 NumPy arrays). This means that when a function that has ``cppmat::cartesian::tensor2s`` as argument, the upper-diagonal part is read; while when it has an argument ``cppmat::cartesian::tensor2d`` only the diagonal is considered.

  **This requires extra attention as information might be lost. To optimize for speed and flexibility no checks are performed in the release libraries derived from cppmat!**

  You can ask cppmat to check for this, by omitting the ``-DNDEBUG`` compiler flag (this enables several assertions, so it may cost you some efficiency).

  (The same holds for the classes and functions from ``<cppmat/tensor2.h>`` and ``<cppmat/tensor3.h>``.)

To use this feature one has to:

.. code-block:: cpp

  #include <cppmat/pybind11.h>

An example is provided in ``docs/examples/tensorlib``. This example includes two forms of building:

1.  ``CMakeList.txt`` for building using ``cmake`` (``cmake .`` and then ``make``). For this to work, *pybind11* must be 'installed' on the system. Alternatively you can include *pybind11* as a sub-folder (for example using ``git submodule add https://github.com/pybind/pybind11.git``). In that case, replace ``find_package(pybind11 REQUIRED)`` by ``add_subdirectory(pybind11)`` in ``CMakeList.txt``.

2.  ``setup.py`` for building using ``python`` (``python setup.py build`` and then ``python setup.py install``). Using this option, ``python`` will take care of the *pybind11* and cppmat dependencies.

    *(Replace the executable with your favorite Python version, e.g. with 'python3')*


