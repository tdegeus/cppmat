
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

  This library is free to use under the `MIT license <https://github.com/tdegeus/cppmat/blob/master/LICENSE>`_. Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word of mouth advertisement, .... Bugs or feature requests can be filed on `GitHub <http://github.com/tdegeus/cppmat>`_. As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

.. tip::

  This document should be considered as a quick-start guide. A lot effort has been spent on the readability of the code itself. One is highly encouraged to answer more advanced questions that arise from this guide directly using the code.

This header-only module provides C++ classes and several accompanying methods to work with n-d arrays and/or tensors. It's usage, programmatically and from a compilation perspective, is really simple. One just has to ``#include <cppmat/cppmat.h>`` and tell your compiler where cppmat is located (and to use the C++14 or younger standard). Really, that's it!

Overview
========

The following dynamically sized classes can be used.

+-------------------------------------------------------------+----------------------------------+
| **Class**                                                   | **Description**                  |
+=============================================================+==================================+
| :ref:`cppmat::array <var_regular_array>`                    | array of arbitrary rank          |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::matrix <var_regular_matrix>`                  | matrix (array of rank 2)         |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::vector <var_regular_vector>`                  | vector (array of rank 1)         |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::symmetric::matrix <var_symmetric_matrix>`     | symmetric, square, matrix        |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::diagonal::matrix <var_diagonal_matrix>`       | diagonal, square, matrix         |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::cartesian::tensor4 <var_cartesian_tensor4>`   | 4th-order tensor                 |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::cartesian::tensor2 <var_cartesian_tensor2>`   | 2nd-order tensor                 |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::cartesian::tensor2s <var_cartesian_tensor2s>` | 2nd-order symmetric tensor       |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::cartesian::tensor2d <var_cartesian_tensor2d>` | 2nd-order diagonal tensor        |
+-------------------------------------------------------------+----------------------------------+
| :ref:`cppmat::cartesian::vector <var_cartesian_vector>`     | 1st-order tensor (a.k.a. vector) |
+-------------------------------------------------------------+----------------------------------+

Each of these classes has a fixed size equivalent (that is usually more efficient):

+-------------------------------------------------------------------+-------------------------------------------------------------+
| **Fixed size**                                                    | **Dynamical size**                                          |
+===================================================================+=============================================================+
| :ref:`cppmat::tiny::array <fix_regular_array>`                    | :ref:`cppmat::array <var_regular_array>`                    |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::matrix <fix_regular_matrix>`                  | :ref:`cppmat::matrix <var_regular_matrix>`                  |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::vector <fix_regular_vector>`                  | :ref:`cppmat::vector <var_regular_vector>`                  |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::symmetric::matrix <fix_symmetric_matrix>`     | :ref:`cppmat::symmetric::matrix <var_symmetric_matrix>`     |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::diagonal::matrix <fix_diagonal_matrix>`       | :ref:`cppmat::diagonal::matrix <var_diagonal_matrix>`       |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::cartesian::tensor4 <fix_cartesian_tensor4>`   | :ref:`cppmat::cartesian::tensor4 <var_cartesian_tensor4>`   |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::cartesian::tensor2 <fix_cartesian_tensor2>`   | :ref:`cppmat::cartesian::tensor2 <var_cartesian_tensor2>`   |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::cartesian::tensor2s <fix_cartesian_tensor2s>` | :ref:`cppmat::cartesian::tensor2s <var_cartesian_tensor2s>` |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::cartesian::tensor2d <fix_cartesian_tensor2d>` | :ref:`cppmat::cartesian::tensor2d <var_cartesian_tensor2d>` |
+-------------------------------------------------------------------+-------------------------------------------------------------+
| :ref:`cppmat::tiny::cartesian::vector <fix_cartesian_vector>`     | :ref:`cppmat::cartesian::vector <var_cartesian_vector>`     |
+-------------------------------------------------------------------+-------------------------------------------------------------+

Each fixed size class has an equivalent which can view a ``const``-pointer (with limited functionality):

+-------------------------------------------------------------------+-------------------------------------------------------------------+
| **View pointer**                                                  | **Fixed size**                                                    |
+===================================================================+===================================================================+
| :ref:`cppmat::view::array <map_regular_array>`                    | :ref:`cppmat::tiny::array <fix_regular_array>`                    |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::matrix <map_regular_matrix>`                  | :ref:`cppmat::tiny::matrix <fix_regular_matrix>`                  |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::vector <map_regular_vector>`                  | :ref:`cppmat::tiny::vector <fix_regular_vector>`                  |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::symmetric::matrix <map_symmetric_matrix>`     | :ref:`cppmat::tiny::symmetric::matrix <fix_symmetric_matrix>`     |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::diagonal::matrix <map_diagonal_matrix>`       | :ref:`cppmat::tiny::diagonal::matrix <fix_diagonal_matrix>`       |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::cartesian::tensor4 <map_cartesian_tensor4>`   | :ref:`cppmat::tiny::cartesian::tensor4 <fix_cartesian_tensor4>`   |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::cartesian::tensor2 <map_cartesian_tensor2>`   | :ref:`cppmat::tiny::cartesian::tensor2 <fix_cartesian_tensor2>`   |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::cartesian::tensor2s <map_cartesian_tensor2s>` | :ref:`cppmat::tiny::cartesian::tensor2s <fix_cartesian_tensor2s>` |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::cartesian::tensor2d <map_cartesian_tensor2d>` | :ref:`cppmat::tiny::cartesian::tensor2d <fix_cartesian_tensor2d>` |
+-------------------------------------------------------------------+-------------------------------------------------------------------+
| :ref:`cppmat::view::cartesian::vector <map_cartesian_vector>`     | :ref:`cppmat::tiny::cartesian::vector <fix_cartesian_vector>`     |
+-------------------------------------------------------------------+-------------------------------------------------------------------+


Example
=======

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

Contents
========

.. toctree::
   :maxdepth: 1

   cppmat_var_regular.rst
   cppmat_var_symmetric.rst
   cppmat_var_diagonal.rst
   cppmat_cartesian.rst
   cppmat_fix.rst
   cppmat_map.rst
   copy.rst
   misc.rst
   histogram.rst
   compile.rst
   python.rst
   develop.rst

