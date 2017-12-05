
.. _cartesian2:

******************
cppmat::cartesian2
******************

.. note::

  See `tensor2.h <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor2.h>`_ and `tensor2.cpp <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor.cpp>`_

A specialization of :ref:`cartesian`, that takes advantage of the knowledge that the arrays are fixed-size and small. Also several loops are unrolled. All the functionality and names of :ref:`cartesian` are transferable to :ref:`cartesian2` and :ref:`cartesian3`.

Compared to :ref:`cartesian` the dimension argument must be omitted everywhere.

.. note::

  There is no way to automatically switch between :ref:`cartesian2`, :ref:`cartesian3`, and :ref:`cartesian`.

Classes
=======

.. _cartesian2_tensor4:

cppmat::cartesian2::tensor4
---------------------------

4th-order tensor (rank 4 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian2::tensor4<double> A;

  A(0,0,0,0) = ...

.. _cartesian2_tensor2:

cppmat::cartesian2::tensor2
---------------------------

2nd-order tensor (rank 2 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian2::tensor2<double> A;

  A(0,0) = ...

.. _cartesian2_tensor2s:

cppmat::cartesian2::tensor2s
----------------------------

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian2::tensor2s<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X , X ;
        X ]

*The remaining components are inferred from symmetry*.

.. _cartesian2_tensor2d:

cppmat::cartesian2::tensor2d
----------------------------

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian2::tensor2d<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X     ;
        X ]

*The remaining components are imposed to be zero*.

.. _cartesian2_vector:

cppmat::cartesian2::vector
--------------------------

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

.. _cartesian3:

******************
cppmat::cartesian3
******************

.. note::

  See `tensor3.h <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor3.h>`_ and `tensor3.cpp <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor3.cpp>`_

Identical to :ref:`cartesian2`, but for 3 dimensions.

Classes
=======

.. _cartesian3_tensor4:

cppmat::cartesian3::tensor4
---------------------------

4th-order tensor (rank 4 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian3::tensor4<double> A;

  A(0,0,0,0) = ...

.. _cartesian3_tensor2:

cppmat::cartesian3::tensor2
---------------------------

2nd-order tensor (rank 2 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian3::tensor2<double> A;

  A(0,0) = ...

.. _cartesian3_tensor2s:

cppmat::cartesian3::tensor2s
----------------------------

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian3::tensor2s<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X , X , X ;
        X , X ;
            X ]

*The remaining components are inferred from symmetry*.

.. _cartesian3_tensor2d:

cppmat::cartesian3::tensor2d
----------------------------

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian3::tensor2d<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X         ;
        X     ;
            X ]

*The remaining components are imposed to be zero*.

.. _cartesian3_vector:

cppmat::cartesian3::vector
--------------------------

Vector (rank 1 tensor) of arbitrary dimension. For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A;

  A(0) = ...
