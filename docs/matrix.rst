
.. _matrix:

*****************
<cppmat/matrix.h>
*****************

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

  If you know that you will work exclusively with a 2-dimensional matrix, please consider using :ref:`matrix2` instead of :ref:`matrix`. This is generally more efficient. If, on top of that, you want to use only a very small matrix, please the fixed sized :ref:`tiny_matrix2`. The latter doesn't need dynamic memory allocation, and can therefore be considerably faster.

Methods
=======

*   ``A(i,j)``

    Return the entry at ``(i,j)``. The number of indices (i.e. ``A(i)``, ``A(i,j)``, ``A(i,j,k)``, ...) has to correspond to the number of dimensions. See :ref:`matrix-index`.

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

    Has the result that ``A.shape() == { 10 , 10 , 1 }``.

.. _matrix-index:

Indexing
========

Periodicity
-----------

Periodic indices can be used. For example ``-1`` will refer to the last index along that dimension (i.e. ``-1 -> N-1``, with ``N = A.shape(i)``), while ``N`` will refer the the first index (``N -> 0``).

This does require a check and possible modification for each index reference. If you know that your indices are within the bounds of the matrix you can avoid it by using indices of type ``size_t``. I.e.

.. code-block:: cpp

  cppmat::matrix<double> A({10,10});

  for ( size_t i = 0 ; i < A.shape(0) ; ++i )
    for ( size_t j = 0 ; j < A.shape(1) ; ++j )
      A(i,j) = ...

is faster than

.. code-block:: cpp

  cppmat::matrix<double> A({10,10});

  for ( int i = 0 ; i < A.shape(0) ; ++i )
    for ( int j = 0 ; j < A.shape(1) ; ++j )
      A(i,j) = ...

since in the latter case the indices are checked, while this check is omitted in the former case.

Number of indices
-----------------

In principle the number of indices should match the dimensions of the matrix (i.e. ``A.ndim()`` and ``A.shape().size()``), though it is no problem to reference to a matrix certain index using a higher-dimensional equivalent. For example:

.. code-block:: cpp

  cppmat::matrix<double> A({10,10});

  A(5,5,0) = ...

is perfectly acceptable. Note that higher-dimensions can only be trailing ones, using for example ``A(0,5,5)`` is not acceptable, nor is of course ``A(5,5,1)``.

View
====

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``
