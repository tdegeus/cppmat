
.. _cppmat:

******
cppmat
******

.. _matrix:

cppmat::matrix
==============

Header-only module that provides a C++ class for n-d matrices. For example a rank 3 matrix is allocated as follows:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

.. note:: **Tip**

  If you know that you will work exclusively with a 2-dimensional matrix, please consider using :ref:`matrix2` instead of :ref:`matrix`. This is generally more efficient. If, on top of that, you want to use only a very small matrix, please the fixed sized :ref:`tiny_matrix2`. The latter doesn't need dynamic memory allocation, and can therefore be considerably faster. For the sake of generality there exist also the classes :ref:`vector` and :ref:`tiny_vector`, but generally one can just just the standard C++ ``std::vector``.

Methods
-------

*   ``A(i,j,k)``

    Return the entry at ``(i,j,k)``. The number of indices (i.e. ``A(i)``, ``A(i,j)``, ``A(i,j,k)``, ...) generally corresponds to the number of dimensions (see :ref:`matrix-index` for additional directives).

*   ``A[i]``

    Returns the ``i``-th entry in the matrix viewed as a list.

*   ``size_t n = A.shape(i)``

    Returns the shape along dimension ``i``.

*   ``std::vector<size_t> shape = A.shape()``

    Returns the shape along all dimensions.

*   ``A.reshape({...})``

    Change the shape of the matrix. It is required that the total number of entries does not change.

*   ``A.chdim(N)``

    Change the number of dimensions to ``N``. This affects the outputted ``shape``. For example:

    .. code-block:: cpp

      cppmat::matrix<double> A({10,10});

      A.chdim(3);

    Has the result that ``A.shape() == {10,10,1}``.

.. _matrix-index:

Indexing
--------

In principle the number of indices should match the dimensions of the matrix (i.e. ``A.ndim()`` and ``A.shape().size()``), though it is no problem to reference to a matrix certain index using a higher-dimensional equivalent. For example:

.. code-block:: cpp

  cppmat::matrix<double> A({10,10});

  A(5,5,0) = ...

is perfectly acceptable. Note that higher-dimensions can only be trailing ones, using for example ``A(0,5,5)`` is not acceptable, nor is of course ``A(5,5,1)``.

Similarly to refer to the beginning of a block (e.g. a row) one can omit the zero arguments. For example to the beginning of the second row of the above matrix one can use ``&A(1)``.

View
----

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``

.. _matrix2:

cppmat::matrix2
===============

Header-only module that provides a C++ class for 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix2<double> A(10,10);

      A(0,0) = ...

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`matrix`, though there is obviously no ``chdim`` method.

.. _vector:

cppmat::vector
==============

Header-only module that provides a C++ class for 1-d matrices (a.k.a. vectors). For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::vector<double> A(10);

      A(0) = ...

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`matrix`, though there is obviously no ``chdim`` method.

.. note::

  Compared to `std::vector` this class is not so much different, with the exception that it provides indexing also with round brackets, and and automated printing of entries.

