
.. _tiny:

************
cppmat::tiny
************

.. _tiny_matrix:

cppmat::tiny::matrix
====================

[:download:`tiny_matrix.h <../src/cppmat/tiny_matrix.h>`, :download:`tiny_matrix.cpp <../src/cppmat/tiny_matrix.cpp>`]

Class for fixed size, small, 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::matrix<double,10,10> A;

      A(0,0) = ...

      ...

      return 0;
  }

Compared to :ref:`regular_matrix` the size of the matrix cannot be dynamically changed. Consequently there is not dynamic memory allocation, often resulting in faster behavior.

.. note::

  The methods are the same as for :ref:`regular_matrix`.

.. note::

  To 'view' a pointer as a matrix, use :ref:`view_tiny_matrix`.

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

  The methods are the same as for :ref:`regular_vector`.

.. note::

  To 'view' a pointer as a matrix, use :ref:`view_tiny_vector`.
