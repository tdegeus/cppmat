
******
cppmat
******

This header-only module provides C++ classes and several accompanying methods to work with n-d matrices and/or tensors. It's usage, programmatically and from a compilation perspective, is really simple. One just has to include the relevant header file and tell your compiler where it is located (and to the C++14 or younger standard). Really, that's it!

After ``#include <cppmat/cppmat.h>`` the following classes can be used.

+--------------------------------------+-------------------------------------------+---------------------+
| **Class**                            | **Description**                           | **Header**          |
+======================================+===========================================+=====================+
| ``cppmat::matrix<T>``                | n-d matrix of flexible size               | :ref:`matrix`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::matrix2<T>``               | 2-d matrix of flexible size               | :ref:`matrix2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::vector<T>``                | array (1-d matrix) of flexible size       | :ref:`vector`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::tiny::matrix2<T,M,N>``     | small, fixed size, 2-d matrix             | :ref:`tiny_matrix2` |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::tiny::vector<T,N>``        | small, fixed size, array                  | :ref:`tiny_vector`  |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian::tensor4<T>``    | n-d tensors of rank 4                     | :ref:`tensor`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian::tensor2<T>``    | n-d tensors of rank 2                     | :ref:`tensor`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian::tensor2s<T>``   | n-d symmetric tensors of rank 2           | :ref:`tensor`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian::tensor2d<T>``   | n-d diagonal tensors of rank 2            | :ref:`tensor`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian::vector<T>``     | n-d diagonal vectors (tensors of rank 1)  | :ref:`tensor`       |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian2d::tensor4<T>``  | 2-d tensors of rank 4                     | :ref:`tensor2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian2d::tensor2<T>``  | 2-d tensors of rank 2                     | :ref:`tensor2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian2d::tensor2s<T>`` | 2-d symmetric tensors of rank 2           | :ref:`tensor2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian2d::tensor2d<T>`` | 2-d diagonal tensors of rank 2            | :ref:`tensor2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian2d::vector<T>``   | 2-d diagonal vectors (tensors of rank 1)  | :ref:`tensor2`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian3d::tensor4<T>``  | 3-d tensors of rank 4                     | :ref:`tensor3`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian3d::tensor2<T>``  | 3-d tensors of rank 2                     | :ref:`tensor3`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian3d::tensor2s<T>`` | 3-d symmetric tensors of rank 2           | :ref:`tensor3`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian3d::tensor2d<T>`` | 3-d diagonal tensors of rank 2            | :ref:`tensor3`      |
+--------------------------------------+-------------------------------------------+---------------------+
| ``cppmat::cartesian3d::vector<T>``   | 3-d diagonal vectors (tensors of rank 1)  | :ref:`tensor3`      |
+--------------------------------------+-------------------------------------------+---------------------+


For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

.. note::

  This library is free to use under the `MIT license <https://github.com/tdegeus/cppmat/blob/master/LICENSE>`_. Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bug reports or feature requests can be filed on `GitHub <http://github.com/tdegeus/cppmat>`_. As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

  Download: `.zip file <https://github.com/tdegeus/cppmat/zipball/master>`_ | `.tar.gz file <https://github.com/tdegeus/cppmat/tarball/master>`_.

  (c - `MIT <https://github.com/tdegeus/cppmat/blob/master/LICENSE>`_) T.W.J. de Geus (Tom) | tom@geus.me | `www.geus.me <http://www.geus.me>`_ | `github.com/tdegeus/cppmat <http://github.com/tdegeus/cppmat>`_

Contents
========

.. toctree::
   :maxdepth: 1

   matrix.rst
   matrix2.rst
   vector.rst
   tiny_matrix2.rst
   tiny_vector.rst
   tensor.rst
   tensor2.rst
   tensor3.rst
   map.rst
   eigen.rst
   compile.rst
   python.rst
   develop.rst

