
.. _cartesian:

*****************
cppmat::cartesian
*****************

.. note::

  See `tensor.h <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor.h>`_ and `tensor.cpp <https://github.com/tdegeus/cppmat/blob/master/src/cppmat/tensor.cpp>`_

Provides classes for 4th- and 2nd order tensors and vectors (the latter essentially coincide with ``std::vector``, but with special methods). For example, a fourth order identity tensor in 3-D is obtained as follows:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    cppmat::cartesian::tensor4<double> A = cppmat::cartesian::identity4<double>(3);

    ...

    return 0;
  }

.. note:: **Tip**

  To print, use the common C++ ``std::cout << A << std::endl;``. To customize formatting use the more classic C syntax ``A.printf("%16.8e");``

.. note:: **Tip**

  If you know that you will work exclusively in 2 or 3 dimensions, please consider using :ref:`cartesian2` or :ref:`cartesian3` instead of :ref:`cartesian`. This is generally more efficient as it can take advantage of the knowledge that the arrays are fixed size and relatively small. Also several loops are unrolled.

.. note:: **Tip**

  The notation can be shortened to:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    using     T4 = cppmat::cartesian::tensor4<double>;
    namespace cm = cppmat::cartesian;

    int main()
    {
      T4 A = cm::identity4<double>(3);

      ...

      return 0;
    }

Classes
=======

.. _cartesian_tensor4:

``cppmat::cartesian::tensor4``
------------------------------

4th-order tensor (rank 4 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian::tensor4<double> A(3);

  A(0,0,0,0) = ...

.. _cartesian_tensor2:

``cppmat::cartesian::tensor2``
------------------------------

2nd-order tensor (rank 2 tensor) of arbitrary dimension.

.. code-block:: cpp

  cppmat::cartesian::tensor2<double> A(3);

  A(0,0) = ...

.. _cartesian_tensor2s:

``cppmat::cartesian::tensor2s``
-------------------------------

Symmetric 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian::tensor2s<double> A(3);

  A(0,0) = ...

For example, for the case of 3 dimensions, the following components are stored:

.. code-block:: cpp

  [ X , X , X ;
        X , X ;
            X ]

*The remaining components are inferred from symmetry*.

.. _cartesian_tensor2d:

``cppmat::cartesian::tensor2d``
-------------------------------

diagonal 2nd-order tensor.

.. code-block:: cpp

  cppmat::cartesian::tensor2d<double> A(3);

  A(0,0) = ...

For example, for the case of 3 dimensions, the following components are stored:

.. code-block:: cpp

  [ X         ;
        X     ;
            X ]

*The remaining components are imposed to be zero*.

.. _cartesian_vector:

``cppmat::cartesian::vector``
-----------------------------

Vector (rank 1 tensor) of arbitrary dimension. For example:

.. code-block:: cpp

  cppmat::cartesian::vector<double> A(3);

  A(0) = ...

.. note::

  Because of the flexibility of C++ it is easy to switch between these specialized classes and the more general ``cppmat::cartesian::tensor2`` classes. For example, the following will work:

  .. code-block:: cpp

    cppmat::cartesian::tensor2d<double> I = cppmat::cartesian::identity2<double>(3);

    cppmat::cartesian::tensor2 <double> A = I;

  or even

  .. code-block:: cpp

    cppmat::cartesian::tensor2 <double> I = cppmat::cartesian::identity2<double>(3);

  Also arithmetic works:

  .. code-block:: cpp

    cppmat::cartesian::tensor2d<double> A = 3.0 * I;

  Note that it is even possible to perform arithmetic between the three different 2nd-order tensor classes, whereby the output type depends on the type of operator.

  Finally, all the :ref:`tensor-methods` accept all three classes - ``cppmat::cartesian::tensor2``, ``cppmat::cartesian::tensor2s``, ``cppmat::cartesian::tensor2d`` - allowing their usage without any prior type casting. In fact the methods will often perform better for the specialized classes since fewer operations are needed.

.. note::

  The easy automatic conversion described above is not possible from a class to another where more assumptions on the structure are made (e.g. from ``cppmat::cartesian::tensor2`` to ``cppmat::cartesian::tensor2d``) because information is (potentially) lost. To still move forward with the conversion, the following manual conversion can be used:

  .. code-block:: cpp

    cppmat::cartesian::tensor2<double> A(3);

    A(0,0) = ...

    // take the symmetric part of "A": "C = (A+A.T())/2."
    cppmat::cartesian::tensor2s<double> C = A.cast<cppmat::cartesian::tensor2s<double>>();

    // take the diagonal of "A"
    cppmat::cartesian::tensor2d<double> C = A.cast<cppmat::cartesian::tensor2d<double>>();

.. note::

  The feature of automatic casting of ``cppmat::cartesian::tensor2s`` and ``cppmat::cartesian::tensor2d`` to a more general class can be switched off by ``#define CPPMAT_NOCONVERT``. It should be remarked that the Python interface relies on this feature.

.. _tensor-methods:

Methods
=======

For each class the index operator ``(...)``, the arithmetic operators ``*=``, ``*``,``/=``, ``/``,``+=``, ``+``,``-=``, ``-``, and the comparison operator ``==`` are available. Also, one can use ``.zeros()`` or ``.ones()`` to initialize all components respectively to zeros or ones. Furthermore, the following methods are available.

.. note::

  Below the rank can be inferred from the indices, but should be easy to understand even without them. Pseudo-code is used to introduce the methods. For the first method it is short for:

  .. code-block:: cpp

    cppmat::cartesian::tensor4<double> A = cppmat::cartesian::identity4<double>(3);
    cppmat::cartesian::tensor2<double> B = cppmat::cartesian::identity2<double>(3);

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

  One can also call the methods as functions using ``cppmmat::ddot( A , B )``, ``cppmmat::dot( A , B )``, ``cppmmat::dyadic( A , B )``, ``cppmmat::cross( A , B )``, ``cppmmat::transpose( A )``, ``cppmmat::transposeR( A )``, ``cppmmat::transposeL( A )``, ``cppmmat::inv( A )``, ``cppmmat::det( A )``, and ``cppmmat::trace( A )``, These methods are just a front end for the class methods described above.

