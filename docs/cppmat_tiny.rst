
.. _tiny:

************
cppmat::tiny
************

.. _tiny_matrix2:

cppmat::tiny::matrix2
=====================

.. note::

  See `tiny_matrix2.h <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tiny_matrix2.h>`_ and `tiny_matrix2.cpp <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tiny_matrix2.cpp>`_

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

Map external pointer
--------------------

This class can be used to 'view' and external pointer. This can be useful to refer to a part of a bigger matrix. For example:

.. code-block:: cpp


  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix2<double> container({100,4,2});

      cppmat::tiny::matrix2<double,4,2> view;

      for ( size_t i = 0 ; i < container.shape(0) ; ++i )
      {
          view.map(&container(i));

          view(0,0) = ... // directly stored in "container"
      }
  }

.. _tiny_vector:

cppmat::tiny::vector
====================

.. note::

  See `tiny_vector.h <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tiny_vector.h>`_ and `tiny_vector.cpp <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tiny_vector.cpp>`_

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
