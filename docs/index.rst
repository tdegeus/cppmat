
******
cppmat
******

.. image:: https://img.shields.io/badge/license-MIT-brightgreen.svg
  :target: https://github.com/tdegeus/cppmat/blob/master/LICENSE
  :alt: MIT license

.. image:: https://img.shields.io/badge/warranty-no-red.svg
  :target: https://github.com/tdegeus/cppmat/blob/master/LICENSE
  :alt: MIT license

.. image:: https://img.shields.io/badge/download-.zip-lightgray.svg
  :target: https://github.com/tdegeus/cppmat/zipball/master
  :alt: Download as .zip

.. image:: https://img.shields.io/badge/download-.tar.gz-lightgray.svg
  :target: https://github.com/tdegeus/cppmat/tarball/master
  :alt: Download as .tar.gz

.. image:: https://img.shields.io/badge/contact-tom@geus.me-blue.svg
  :target: mailto:tom@geus.me
  :alt: Contact tom@geus.me

.. image:: https://img.shields.io/badge/contact-cppmat.geus.me-blue.svg
  :target: http://cppmat.geus.me
  :alt: Website cppmat.geus.me

.. image:: https://img.shields.io/badge/contact-www.geus.me-blue.svg
  :target: http://www.geus.me
  :alt: Website www.geus.me

.. image:: https://img.shields.io/badge/GitHub-tdegeus/cppmat-blue.svg
  :target: https://github.com/tdegeus/cppmat
  :alt: Github tdegeus/cppmat


.. note::

  This document should be considered as a quick-start guide. A lot effort has been spent on the readability of the code itself. One is highly encouraged to answer more advanced questions that arise from this guide directly using the code.

This header-only module provides C++ classes and several accompanying methods to work with n-d matrices and/or tensors. It's usage, programmatically and from a compilation perspective, is really simple. One just has to ``#include <cppmat/cppmat.h>`` and tell your compiler where cppmat is located (and to the C++14 or younger standard). Really, that's it!


The following classes can be used.

+--------------------------------------+-------------------------------------------+----------------------------+
| **Class**                            | **Description**                           | **Header**                 |
+======================================+===========================================+============================+
| ``cppmat::matrix<T>``                | n-d matrix of flexible size               | :ref:`matrix`              |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::matrix2<T>``               | 2-d matrix of flexible size               | :ref:`matrix2`             |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::vector<T>``                | array (1-d matrix) of flexible size       | :ref:`vector`              |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::periodic::matrix<T>``      | n-d matrix of flexible size               | :ref:`periodic_matrix`     |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::periodic::matrix2<T>``     | 2-d matrix of flexible size               | :ref:`periodic_matrix2`    |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::periodic::vector<T>``      | array (1-d matrix) of flexible size       | :ref:`periodic_vector`     |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::tiny::matrix2<T,M,N>``     | small, fixed size, 2-d matrix             | :ref:`tiny_matrix2`        |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::tiny::vector<T,N>``        | small, fixed size, array                  | :ref:`tiny_vector`         |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian::tensor4<T>``    | n-d tensors of rank 4                     | :ref:`cartesian_tensor4`   |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian::tensor2<T>``    | n-d tensors of rank 2                     | :ref:`cartesian_tensor2`   |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian::tensor2s<T>``   | n-d symmetric tensors of rank 2           | :ref:`cartesian_tensor2s`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian::tensor2d<T>``   | n-d diagonal tensors of rank 2            | :ref:`cartesian_tensor2d`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian::vector<T>``     | n-d diagonal vectors (tensors of rank 1)  | :ref:`cartesian_vector`    |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian2d::tensor4<T>``  | 2-d tensors of rank 4                     | :ref:`cartesian2_tensor4`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian2d::tensor2<T>``  | 2-d tensors of rank 2                     | :ref:`cartesian2_tensor2`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian2d::tensor2s<T>`` | 2-d symmetric tensors of rank 2           | :ref:`cartesian2_tensor2s` |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian2d::tensor2d<T>`` | 2-d diagonal tensors of rank 2            | :ref:`cartesian2_tensor2d` |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian2d::vector<T>``   | 2-d diagonal vectors (tensors of rank 1)  | :ref:`cartesian2_vector`   |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian3d::tensor4<T>``  | 3-d tensors of rank 4                     | :ref:`cartesian3_tensor4`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian3d::tensor2<T>``  | 3-d tensors of rank 2                     | :ref:`cartesian3_tensor2`  |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian3d::tensor2s<T>`` | 3-d symmetric tensors of rank 2           | :ref:`cartesian3_tensor2s` |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian3d::tensor2d<T>`` | 3-d diagonal tensors of rank 2            | :ref:`cartesian3_tensor2d` |
+--------------------------------------+-------------------------------------------+----------------------------+
| ``cppmat::cartesian3d::vector<T>``   | 3-d diagonal vectors (tensors of rank 1)  | :ref:`cartesian3_vector`   |
+--------------------------------------+-------------------------------------------+----------------------------+


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

  (c - `MIT <https://github.com/tdegeus/cppmat/blob/master/LICENSE>`_) T.W.J. de Geus (Tom) | tom@geus.me | `www.geus.me <http://www.geus.me>`_ | `github.com/tdegeus/cppmat <http://github.com/tdegeus/cppmat>`_

Contents
========

.. toctree::
   :maxdepth: 1

   cppmat.rst
   cppmat_periodic.rst
   cppmat_tiny.rst
   cppmat_tensor.rst
   cppmat_tensor2.rst
   copy.rst
   compile.rst
   python.rst
   develop.rst

