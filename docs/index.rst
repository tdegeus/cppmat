
******
cppmat
******

.. |badge1| image:: https://img.shields.io/badge/license-MIT-brightgreen.svg
  :target: https://github.com/tdegeus/cppmat/blob/master/LICENSE
  :alt: MIT license

.. |badge2| image:: https://img.shields.io/badge/warranty-no-red.svg
  :target: https://github.com/tdegeus/cppmat/blob/master/LICENSE
  :alt: MIT license

.. |badge3| image:: https://img.shields.io/badge/download-.zip-lightgray.svg
  :target: https://github.com/tdegeus/cppmat/zipball/master
  :alt: Download as .zip

.. |badge4| image:: https://img.shields.io/badge/download-.tar.gz-lightgray.svg
  :target: https://github.com/tdegeus/cppmat/tarball/master
  :alt: Download as .tar.gz

.. |badge5| image:: https://img.shields.io/badge/contact-tom@geus.me-blue.svg
  :target: mailto:tom@geus.me
  :alt: Contact tom@geus.me

.. |badge6| image:: https://img.shields.io/badge/contact-www.geus.me-blue.svg
  :target: http://www.geus.me
  :alt: Website www.geus.me

.. |badge7| image:: https://img.shields.io/badge/GitHub-tdegeus/cppmat-blue.svg
  :target: https://github.com/tdegeus/cppmat
  :alt: Github tdegeus/cppmat

.. |badge8| image:: https://img.shields.io/badge/documentation-cppmat.geus.me-blue.svg
  :target: http://cppmat.geus.me
  :alt: Website cppmat.geus.me

| |badge1| |badge2| |badge3| |badge4|
| |badge5| |badge6| |badge7|
| |badge8|

.. note::

  This library is free to use under the `MIT license <https://github.com/tdegeus/cppmat/blob/master/LICENSE>`_. Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bug reports or feature requests can be filed on `GitHub <http://github.com/tdegeus/cppmat>`_. As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

.. note::

  This document should be considered as a quick-start guide. A lot effort has been spent on the readability of the code itself. One is highly encouraged to answer more advanced questions that arise from this guide directly using the code.

This header-only module provides C++ classes and several accompanying methods to work with n-d matrices and/or tensors. It's usage, programmatically and from a compilation perspective, is really simple. One just has to ``#include <cppmat/cppmat.h>`` and tell your compiler where cppmat is located (and to the C++14 or younger standard). Really, that's it!


Overview
========

The following classes can be used.

+-------------------------------------------------------------+--------------------------------------+
| **Class**                                                   | **Description**                      |
+=============================================================+======================================+
| :ref:`cppmat::matrix <matrix>`                              | n-d matrix of flexible size          |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::matrix2 <matrix2>`                            | 2-d matrix of flexible size          |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::vector <vector>`                              | array (1-d matrix) of flexible size  |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::periodic::matrix <periodic_matrix>`           | n-d matrix of flexible size          |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::periodic::matrix2 <periodic_matrix2>`         | 2-d matrix of flexible size          |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::periodic::vector <periodic_vector>`           | array (1-d matrix) of flexible size  |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::tiny::matrix2 <tiny_matrix2>`                 | small, fixed size, 2-d matrix        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::tiny::vector <tiny_vector>`                   | small, fixed size, array             |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian::tensor4 <cartesian_tensor4>`       | n-d tensor of rank 4                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian::tensor2 <cartesian_tensor2>`       | n-d tensor of rank 2                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian::tensor2s <cartesian_tensor2s>`     | n-d symmetric tensor of rank 2       |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian::tensor2d <cartesian_tensor2d>`     | n-d diagonal tensor of rank 2        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian::vector <cartesian_vector>`         | n-d vector (tensor of rank 1)        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian2d::tensor4 <cartesian2d_tensor4>`   | 2-d tensor of rank 4                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian2d::tensor2 <cartesian2d_tensor2>`   | 2-d tensor of rank 2                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian2d::tensor2s <cartesian2d_tensor2s>` | 2-d symmetric tensor of rank 2       |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian2d::tensor2d <cartesian2d_tensor2d>` | 2-d diagonal tensor of rank 2        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian2d::vector <cartesian2d_vector>`     | 2-d vector (tensor of rank 1)        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian3d::tensor4 <cartesian3d_tensor4>`   | 3-d tensor of rank 4                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian3d::tensor2 <cartesian3d_tensor2>`   | 3-d tensor of rank 2                 |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian3d::tensor2s <cartesian3d_tensor2s>` | 3-d symmetric tensor of rank 2       |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian3d::tensor2d <cartesian3d_tensor2d>` | 3-d diagonal tensor of rank 2        |
+-------------------------------------------------------------+--------------------------------------+
| :ref:`cppmat::cartesian3d::vector <cartesian3d_vector>`     | 3-d vector (tensor of rank 1)        |
+-------------------------------------------------------------+--------------------------------------+

In addition, the following classes are available to view a ``const``-pointer:

+------------------------------------------------------------------------+-------------------------------------------------------------+
| **Class**                                                              | **Equivalent class**                                        |
+========================================================================+=============================================================+
| :ref:`cppmat::view::matrix2 <view_matrix2>`                            | :ref:`cppmat::tiny::matrix2 <tiny_matrix2>`                 |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::vector <view_vector>`                              | :ref:`cppmat::tiny::vector <tiny_vector>`                   |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian2d::tensor4 <view_cartesian2d_tensor4>`   | :ref:`cppmat::cartesian2d::tensor4 <cartesian2d_tensor4>`   |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian2d::tensor2 <view_cartesian2d_tensor2>`   | :ref:`cppmat::cartesian2d::tensor2 <cartesian2d_tensor2>`   |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian2d::tensor2s <view_cartesian2d_tensor2s>` | :ref:`cppmat::cartesian2d::tensor2s <cartesian2d_tensor2s>` |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian2d::tensor2d <view_cartesian2d_tensor2d>` | :ref:`cppmat::cartesian2d::tensor2d <cartesian2d_tensor2d>` |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian2d::vector <view_cartesian2d_vector>`     | :ref:`cppmat::cartesian2d::vector <cartesian2d_vector>`     |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian3d::tensor4 <view_cartesian3d_tensor4>`   | :ref:`cppmat::cartesian3d::tensor4 <cartesian3d_tensor4>`   |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian3d::tensor2 <view_cartesian3d_tensor2>`   | :ref:`cppmat::cartesian3d::tensor2 <cartesian3d_tensor2>`   |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian3d::tensor2s <view_cartesian3d_tensor2s>` | :ref:`cppmat::cartesian3d::tensor2s <cartesian3d_tensor2s>` |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian3d::tensor2d <view_cartesian3d_tensor2d>` | :ref:`cppmat::cartesian3d::tensor2d <cartesian3d_tensor2d>` |
+------------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::view::cartesian3d::vector <view_cartesian3d_vector>`     | :ref:`cppmat::cartesian3d::vector <cartesian3d_vector>`     |
+------------------------------------------------------------------------+-------------------------------------------------------------+

Example
=======

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

Contents
========

.. toctree::
   :maxdepth: 1

   cppmat.rst
   cppmat_periodic.rst
   cppmat_tiny.rst
   cppmat_view.rst
   cppmat_cartesian.rst
   cppmat_cartesian2d.rst
   cppmat_cartesian3d.rst
   cppmat_view_cartesian2d.rst
   cppmat_view_cartesian3d.rst
   copy.rst
   compile.rst
   python.rst
   develop.rst

