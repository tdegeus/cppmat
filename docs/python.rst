
.. _python:

****************
Python interface
****************

This library includes provides an interface to `pybind11 <https://github.com/pybind/pybind11>`_ such that an interface to *NumPy* arrays is automatically provided when including a function with a ``cppmat::matrix`` (rank n *NumPy*-array), ``cppmat::tensor4`` (rank 4 *NumPy*-array), ``cppmat::tensor2`` (rank 2 *NumPy*-array), ``cppmat::tensor2s`` (rank 2 *NumPy*-array), ``cppmat::tensor2d`` (rank 2 *NumPy*-array), or ``cppmat::vector`` (rank 1 *NumPy*-array). To use this feature one has to include (either or both):

.. code-block:: cpp

  #include <cppmat/pybind11_matrix.h>
  #include <cppmat/pybind11_tensor.h>

An example is provided in ``docs/examples/tensorlib``. This example includes two forms of building:

1.  ``CMakeList.txt`` for building using ``cmake`` (``cmake .`` and then ``make``). For this to work, *pybind11* must be 'installed' on the system. Alternatively you can include *pybind11* as a sub-folder (for example using ``git submodule add https://github.com/pybind/pybind11.git``). In that case, replace ``find_package(pybind11 REQUIRED)`` by ``add_subdirectory(pybind11)`` in ``CMakeList.txt``.

2.  ``setup.py`` for building using ``python`` (`python setup.py build` and then `python setup.py install`). Using this option, ``python`` will take care of the *pybind11* and *cppmat* dependencies.

    *(Replace the executable with your favorite Python version, e.g. with ``python3``)*

.. note:: **Warning**

  On the Python side all 2nd-order tensors (``cppmat::tensor2``, ``cppmat::tensor2s``, and ``cppmat::tensor2d``) are the same rank 2 *NumPy*-array. This means that when a function with has ``cppmat::tensor2s`` as argument, the upper-diagonal part is read; while when it has an argument ``cppmat::tensor2d`` only the diagonal is considered. You can ask *cppmat* to check for this, by omitting the ``-DNDEBUG`` compiler flag (this enables several assertions, so it may cost you some efficiency).

  **This requires extra attention as information might be lost. To optimize for speed and flexibility no checks are performed in the release libraries derived from *cppmat*!**

