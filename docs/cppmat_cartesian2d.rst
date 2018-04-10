
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

4th-order tensor (rank 4 tensor).

.. code-block:: cpp

  cppmat::cartesian2d::tensor4<double> A;

  A(0,0,0,0) = ...

.. _cartesian2d_tensor2:

cppmat::cartesian2d::tensor2
----------------------------

2nd-order tensor (rank 2 tensor).

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

Diagonal 2nd-order tensor.

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

Vector (rank 1 tensor). For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A;

  A(0) = ...
