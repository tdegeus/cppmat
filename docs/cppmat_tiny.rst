
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

Compared to :ref:`matrix2`:

*   The size of the matrix cannot be dynamically changed.

*   This class can be used to 'view' and external pointer.

The rest of the interface is the same.

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

Compared to :ref:`vector`:

*   The size of the matrix cannot be dynamically changed.

*   This class can be used to 'view' and external pointer.

The rest of the interface is the same.

.. note::

  To 'view' a pointer as a vector, use :ref:`view_vector`.

