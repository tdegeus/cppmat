
.. _cartesian2d:

*******************
cppmat::cartesian2d
*******************

[:download:`tensor2.h <../src/cppmat/tensor2.h>`, :download:`tensor2.cpp <../src/cppmat/tensor2.cpp>`]

A specialization of :ref:`cartesian`, that takes advantage of the knowledge that the arrays are fixed size and small. Also several loops are unrolled. All the functionality and names of :ref:`cartesian` are transferable to :ref:`cartesian2d` and :ref:`cartesian3d`.

Compared to :ref:`cartesian` the dimension argument must be omitted everywhere.

.. note::

  There is no way to automatically switch between :ref:`cartesian2d`, :ref:`cartesian3d`, and :ref:`cartesian`.

Classes
=======

.. _cartesian2d_tensor4:

cppmat::cartesian2d::tensor4
----------------------------

4th-order tensor (rank 4 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian2d::tensor4<double> A;

  A(0,0,0,0) = ...

.. _cartesian2d_tensor2:

cppmat::cartesian2d::tensor2
----------------------------

2nd-order tensor (rank 2 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian2d::tensor2<double> A;

  A(0,0) = ...

.. _cartesian2d_tensor2s:

cppmat::cartesian2d::tensor2s
-----------------------------

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian2d::tensor2s<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X , X ;
        X ]

*The remaining components are inferred from symmetry*.

.. _cartesian2d_tensor2d:

cppmat::cartesian2d::tensor2d
-----------------------------

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian2d::tensor2d<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X     ;
        X ]

*The remaining components are imposed to be zero*.

.. _cartesian2d_vector:

cppmat::cartesian2d::vector
---------------------------

Vector (rank 1 tensor) of arbitrary dimension. For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A;

  A(0) = ...

Map external pointer
====================

Like in :ref:`tiny`, the classes under :ref:`cartesian2` can be used to 'view' an external pointer. For example, for a matrix of 2-d symmetric tensors:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix2<double> container({50,50,3});

      cppmat::cartesian2d::tensor2s view;

      for ( size_t i = 0 ; i < container.shape(0) ; ++i )
      {
          for ( size_t j = 0 ; j < container.shape(1) ; ++j )
          {
              view.map(&container(i,j));

              view(0,0) = ... // directly stored in "container"
          }
      }
  }

.. note::

  The situation can occur that you want to map a ``const`` pointer, for example when you are designing a function that reads from a matrix that is marked ``const`` (because you want to use the matrix 'read-only'). In that case the ``cppmat::cartesian2d`` tensor that you use to map the larger object should be templated using ``const`` (e.g. ``const double``, ``const size_t``, etc.). The ``cppmat::cartesian2d`` tensor can then only be used 'read-only'.

  This case is illustrated in this example:

  .. code-block:: cpp

    #include "cppmat.h"

    void view(const cppmat::matrix<double> &matrix)
    {
      cppmat::cartesian2d::tensor2s<const double> view;

      for ( auto i = 0 ; i < matrix.shape(0) ; ++i )
      {
        view.map(&matrix(i));
        std::cout << view << std::endl;
      }
    }


    int main()
    {
      cppmat::matrix<double> A({2,3});

      A(0,0) = 1.0;  A(1,0) = 101.0;
      A(0,1) = 2.0;  A(1,1) = 102.0;
      A(0,2) = 3.0;  A(1,2) = 103.0;

      view(A);

      return 0;
    }
