
.. _var_cartesian:

*****************
cppmat::cartesian
*****************

Provides classes for 4th- and 2nd order tensors and vectors. For example, a fourth order identity tensor in 3-D is obtained as follows:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    cppmat::cartesian::tensor4<double> A = cppmat::cartesian::tensor4<double>::I(3);

    ...

    std::cout << A << std::endl;

    return 0;
  }

.. tip::

  If you know that you will work in a fixed (and small) number of dimensions (e.g. 2 or 3), please consider using :ref:`fix_cartesian` instead of :ref:`var_cartesian`. This is generally more efficient as it can take advantage of the knowledge that the arrays are fixed size and relatively small. Also several loops are unrolled.

.. tip::

  The notation can be shortened to:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    using T4 = cppmat::cartesian::tensor4<double>;

    int main()
    {
      T4 A = T4::I(3);

      ...

      return 0;
    }

Classes
=======

.. _var_cartesian_tensor4:

``cppmat::cartesian::tensor4``
------------------------------

[:download:`var_cartesian_tensor4.h <../src/cppmat/var_cartesian_tensor4.h>`, :download:`var_cartesian_tensor4.hpp <../src/cppmat/var_cartesian_tensor4.hpp>`]

4th-order tensor (rank 4 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian::tensor4<double> A(3);

  A(0,0,0,0) = ...

.. _var_cartesian_tensor2:

``cppmat::cartesian::tensor2``
------------------------------

[:download:`var_cartesian_tensor2.h <../src/cppmat/var_cartesian_tensor2.h>`, :download:`var_cartesian_tensor2.hpp <../src/cppmat/var_cartesian_tensor2.hpp>`]

2nd-order tensor (rank 2 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian::tensor2<double> A(3);

  A(0,0) = ...

.. _var_cartesian_tensor2s:

``cppmat::cartesian::tensor2s``
-------------------------------

[:download:`var_cartesian_tensor2s.h <../src/cppmat/var_cartesian_tensor2s.h>`, :download:`var_cartesian_tensor2s.hpp <../src/cppmat/var_cartesian_tensor2s.hpp>`]

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian::tensor2s<double> A(3);

  A(0,0) = ...

For example, for the case of 3 dimensions, the following components are stored:

.. code-block:: cpp

  [ X , X , X ;
        X , X ;
            X ]

*The remaining components are inferred from symmetry*. See :ref:`var_symmetric_matrix`.

.. _var_cartesian_tensor2d:

``cppmat::cartesian::tensor2d``
-------------------------------

[:download:`var_cartesian_tensor2d.h <../src/cppmat/var_cartesian_tensor2d.h>`, :download:`var_cartesian_tensor2d.hpp <../src/cppmat/var_cartesian_tensor2d.hpp>`]

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian::tensor2d<double> A(3);

  A(0,0) = ...

For example, for the case of 3 dimensions, the following components are stored:

.. code-block:: cpp

  [ X         ;
        X     ;
            X ]

*The remaining components are imposed to be zero*. See :ref:`var_diagonal_matrix`.

.. _var_cartesian_vector:

``cppmat::cartesian::vector``
-----------------------------

[:download:`var_cartesian_vector.h <../src/cppmat/var_cartesian_vector.h>`, :download:`var_cartesian_vector.hpp <../src/cppmat/var_cartesian_vector.hpp>`]

Vector (rank 1 tensor) of arbitrary dimension. For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A(3);

  A(0) = ...

.. note::

  Because of the flexibility of C++ it is easy to switch between these specialized classes and the more general ``cppmat::cartesian::tensor2`` classes. For example, the following will work:

  .. code-block:: cpp

    using T2  = cppmat::cartesian::tensor2 <double>;
    using T2d = cppmat::cartesian::tensor2d<double>;

    T2d I = T2d::I(3);
    T2  A = I;

  or even

  .. code-block:: cpp

    T2 I = T2d::I(3);

  Also arithmetic works:

  .. code-block:: cpp

    T2d A = 3.0 * I;

  Note that it is even possible to perform arithmetic between the three different 2nd-order tensor classes, whereby the output type depends on the type of operator.

  Finally, all the :ref:`tensor-methods` accept all three classes - ``cppmat::cartesian::tensor2``, ``cppmat::cartesian::tensor2s``, ``cppmat::cartesian::tensor2d`` - allowing their usage without any prior type casting. In fact the methods will often perform better for the specialized classes since fewer operations are needed.

.. note::

  The easy automatic conversion described above is not possible from a class to another where more assumptions on the structure are made (e.g. from ``cppmat::cartesian::tensor2`` to ``cppmat::cartesian::tensor2d``) because information is (potentially) lost.

.. _tensor-methods:

Methods
=======

All the methods of :ref:`var_regular_array` (or :ref:`var_regular_matrix`, :ref:`var_symmetric_matrix`, :ref:`var_symmetric_matrix`, or :ref:`var_regular_vector`) are overloaded. In addition, the following tensor algebra is available.

.. note::

  Below the rank can be inferred from the indices, but should be easy to understand even without them. Pseudo-code is used to introduce the methods. For the first method it is short for:

  .. code-block:: cpp

    cppmat::cartesian::tensor4<double> A = cppmat::cartesian::tensor4<double>::I(3);
    cppmat::cartesian::tensor2<double> B = cppmat::cartesian::tensor2<double>::I(3);

    cppmat::cartesian::tensor2<double> C = A.ddot(B);

  Finally, each occurrence of ``cppmat::cartesian::tensor2`` can be replaced by ``cppmat::cartesian::tensor2s`` or ``cppmat::cartesian::tensor2d``. The latter two often perform better.

*   ``cppmat::cartesian::tensor4<X>``:

    -   ``cppmat::cartesian::tensor4<X> C = A.ddot(const cppmat::cartesian::tensor4<X> &B)``

        Double tensor contraction : :math:`C_{ijmn} = A_{ijkl} B_{lkmn}`

    -   ``cppmat::cartesian::tensor2<X> C = A.ddot(const cppmat::cartesian::tensor2<X> &B)``

        Double tensor contraction :math:`C_{ij} = A_{ijkl} B_{lk}`

    -   ``cppmat::cartesian::tensor4<X> C = A.T()``

        Transposition :math:`C_{lkji} = A_{ijkl}`

    -   ``cppmat::cartesian::tensor4<X> C = A.LT()``

        Left transposition :math:`C_{jikl} = A_{ijkl}`

    -   ``cppmat::cartesian::tensor4<X> C = A.RT()``

        Right transposition :math:`C_{ijlk} = A_{ijkl}`

*   ``cppmat::cartesian::tensor2<X>``:

    -   ``cppmat::cartesian::tensor2<X> C = A.ddot(const cppmat::cartesian::tensor4<X> &B)``

        Double tensor contraction :math:`C_{kl} = A_{ij} B_{jikl}`

    -   ``X C = A.ddot(const cppmat::cartesian::tensor2<X> &B)``

        Double tensor contraction :math:`C = A_{ij} B_{ji}`

    -   ``cppmat::cartesian::tensor2<X> C = A.dot(const cppmat::cartesian::tensor2<X> &B)``

        Tensor contraction :math:`C_{ik} = A_{ij} B_{jk}`

    -   ``cppmat::cartesian::vector<X> C = A.dot(const cppmat::cartesian::vector<X> &B)``

        Tensor contraction :math:`C_{i} = A_{ij} B_{j}`

    -   ``cppmat::cartesian::tensor4<X> C = A.dyadic(const cppmat::cartesian::tensor2<X> &B)``

        Dyadic product :math:`C_{ijkl} = A_{ij} B_{kl}`

    -   ``cppmat::cartesian::tensor2<X> C = A.T()``

        Transposition :math:`C_{ji} = A_{ij}`

    -   ``X C = A.trace()``

        The trace of the tensor (i.e. the sum of the diagonal components) :math:`C = A_{ii}`

    -   ``X C = A.det()``

        The determinant :math:`C = \det \underline{\bm{A}}`

    -   ``cppmat::cartesian::tensor2<X> C = A.inv()``

        The inverse :math:`C_{ij} = A_{ij}^{-1}`

*   ``cppmat::cartesian::vector<X>``:

    -   ``X C = A.dot(const cppmat::cartesian::vector<X> &B)``

        Tensor contraction :math:`C = A_{i} B_{i}`

    -   ``cppmat::cartesian::vector<X> C = A.dot(const cppmat::cartesian::tensor2<X> &B)``

        Tensor contraction :math:`C_{j} = A_{i} B_{ij}`

    -   ``cppmat::cartesian::tensor2<X> C = A.dyadic(const cppmat::cartesian::vector<X> &B)``

        Dyadic product :math:`C_{ij} = A_{i} B_{j}`

    -   ``cppmat::cartesian::vector<X> C = A.cross(const cppmat::cartesian::vector<X> &B)``

        Cross product :math:`\vec{C} = \vec{A} \otimes \vec{B}`


.. note::

  One can also call the methods as functions using ``cppmmat::ddot(A,B)``, ``cppmmat::dot(A,B)``, ``cppmmat::dyadic(A,B)``, ``cppmmat::cross(A,B)``, ``cppmmat::T(A)``, ``cppmmat::RT(A)``, ``cppmmat::LT(A)``, ``cppmmat::inv(A)``, ``cppmmat::det(A)``, and ``cppmmat::trace(A)``. This is fully equivalent (in fact the class methods call these external functions).

