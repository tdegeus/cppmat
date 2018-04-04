
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

  If you wish to copy from a pointer, rather than to point to an external object, the `.copy(...)` function is available with an identical syntax to `.map(...)`.

.. note::

  To map a ``const``-pointer (read-only):

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
        cppmat::matrix2<double> container({50,50,3});

        cppmat::view::cartesian2d::tensor2s view;

        for ( size_t i = 0 ; i < container.shape(0) ; ++i )
        {
            for ( size_t j = 0 ; j < container.shape(1) ; ++j )
            {
                view.map(&container(i,j));

                std::cout << view << std::endl;
            }
        }
    }
