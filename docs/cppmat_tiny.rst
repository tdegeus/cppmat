
.. _tiny:

************
cppmat::tiny
************

.. _tiny_matrix2:

cppmat::tiny::matrix2
=====================

[:download:`tiny_matrix2.h <../src/cppmat/tiny_matrix2.h>`, :download:`tiny_matrix2.cpp <../src/cppmat/tiny_matrix2.cpp>`]

Class for fixed size, small, 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::matrix2<double,10,10> A;

      A(0,0) = ...

      ...

      return 0;
  }

Compared to :ref:`matrix2` the size of the matrix cannot be dynamically changed. Consequently there is not dynamic memory allocation, often resulting in faster behavior.

.. note::

  The methods are the same as for :ref:`matrix2`.

.. note::

  To 'view' a pointer as a matrix, use :ref:`view_matrix2`.

.. _tiny_vector:

cppmat::tiny::vector
====================

[:download:`tiny_vector.h <../src/cppmat/tiny_vector.h>`, :download:`tiny_vector.cpp <../src/cppmat/tiny_vector.cpp>`]

Class for fixed size, small, 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::vector<double,10> A;

      A(0) = ...

      ...

      return 0;
  }

.. note::

  The methods are the same as for :ref:`vector`.

.. note::

  To 'view' a pointer as a matrix, use :ref:`view_vector`.
