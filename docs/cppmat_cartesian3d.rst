
.. _cartesian3d:

*******************
cppmat::cartesian3d
*******************

[:download:`tensor3.h <../src/cppmat/tensor3.h>`, :download:`tensor3.cpp <../src/cppmat/tensor3.cpp>`]

Identical to :ref:`cartesian2d`, but for 3 dimensions.

Classes
=======

.. _cartesian3d_tensor4:

cppmat::cartesian3d::tensor4
----------------------------

4th-order tensor (rank 4 tensor).

.. code-block:: cpp

  cppmat::cartesian3d::tensor4<double> A;

  A(0,0,0,0) = ...

.. _cartesian3d_tensor2:

cppmat::cartesian3d::tensor2
----------------------------

2nd-order tensor (rank 2 tensor).

.. code-block:: cpp

  cppmat::cartesian3d::tensor2<double> A;

  A(0,0) = ...

.. _cartesian3d_tensor2s:

cppmat::cartesian3d::tensor2s
-----------------------------

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian3d::tensor2s<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X , X , X ;
        X , X ;
            X ]

*The remaining components are inferred from symmetry*.

.. _cartesian3d_tensor2d:

cppmat::cartesian3d::tensor2d
-----------------------------

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian3d::tensor2d<double> A;

  A(0,0) = ...

The following components are stored:

.. code-block:: cpp

  [ X         ;
        X     ;
            X ]

*The remaining components are imposed to be zero*.

.. _cartesian3d_vector:

cppmat::cartesian3d::vector
---------------------------

Vector (rank 1 tensor). For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A;

  A(0) = ...
