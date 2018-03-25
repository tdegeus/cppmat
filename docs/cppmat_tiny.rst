
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

.. note::

  The situation can occur that you want to map a ``const`` pointer, for example when you are designing a function that reads from a matrix that is marked ``const`` (because you want to use the matrix 'read-only'). In that case the ``cppmat::tiny`` matrix (or vector) that you use to map the larger object should be templated using ``const`` (e.g. ``const double``, ``const size_t``, etc.). The ``cppmat::tiny`` matrix can then only be used 'read-only'.

  This case is illustrated in this example:

  .. code-block:: cpp

    #include "cppmat.h"

    void view(const cppmat::matrix<size_t> &matrix)
    {
      cppmat::tiny::matrix2<const size_t,2,2> view;

      for ( auto i = 0 ; i < matrix.shape(0) ; ++i )
      {
        view.map(&matrix(i));
        std::cout << view << std::endl;
      }
    }


    int main()
    {
      cppmat::matrix<size_t> A({2,2,2});

      A(0,0,0) =  1;  A(1,0,0) = 101;
      A(0,0,1) =  2;  A(1,0,1) = 102;
      A(0,1,0) = 11;  A(1,1,0) = 111;
      A(0,1,1) = 12;  A(1,1,1) = 112;

      view(A);

      return 0;
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
