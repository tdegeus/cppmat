
.. _cartesian3:

******************
cppmat::cartesian3
******************

[:download:`tensor3.h <../src/cppmat/tensor3.h>`, :download:`tensor3.cpp <../src/cppmat/tensor3.cpp>`]

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
