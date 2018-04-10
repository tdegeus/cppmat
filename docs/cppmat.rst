
.. _cppmat:

******
cppmat
******

.. _matrix:

cppmat::matrix
==============

[:download:`matrix.h <../src/cppmat/matrix.h>`, :download:`matrix.cpp <../src/cppmat/matrix.cpp>`]

Header-only module that provides a C++ class for n-d matrices. For example, a rank 3 matrix is allocated as follows:

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

  If you know that you will work exclusively with a 2-dimensional matrix, please consider using :ref:`matrix2` instead of :ref:`matrix`. This is generally more efficient. If, on top of that, you want to use only a small matrix, please the fixed sized :ref:`tiny_matrix2`. It doesn't need dynamic memory allocation, and can therefore be considerably faster. For the sake of generality there exist also the classes :ref:`vector` and :ref:`tiny_vector`, but generally one can just use the standard C++ ``std::vector``.

Methods
-------

*   ``A(i,j,k)``

    Return the entry at ``(i,j,k)``. The number of indices (i.e. ``A(i)``, ``A(i,j)``, ``A(i,j,k)``, ...) generally corresponds to the number of dimensions (see :ref:`matrix-index` for additional directives).

*   ``A[i]``

    Returns the ``i``-th entry in the matrix viewed as a list.

*   ``A.at(first, last)``

    Returns the entry ``{i,j,k}``, which are stored in a list. The function takes an iterator to the first and the last index of this list. See :ref:`matrix-index-advanced`.

*   ``A.item(i,j,k)``

    Return an iterator to the entry at ``(i,j,k)``.

*   ``A.begin()``

    Return an iterator to the first entry of the matrix.

*   ``A.end()``

    Return an iterator to the last entry of the matrix.

*   ``size_t n = A.shape(i)``

    Returns the shape along dimension ``i``.

*   ``std::vector<size_t> shape = A.shape()``

    Returns the shape along all dimensions.

*   ``A.resize({...})``

    Resize the matrix.

*   ``A.reshape({...})``

    Change the shape of the matrix. It is required that the total number of entries does not change.

*   ``A.chdim(N)``

    Change the number of dimensions to ``N``. This affects the outputted ``shape``. For example:

    .. code-block:: cpp

      cppmat::matrix<double> A({10,10});

      A.chdim(3);

    Has the result that ``A.shape() == {10,10,1}``.

*   ``A.setConstant(...)``, ``A.setZero(...)``, ``A.setOnes(...)``, ``A.zeros(...)``, ``A.ones(...)``

    Set all entries to a constant, zero, or one.

.. _matrix-index:

Indexing
--------

In principle the number of indices should match the dimensions of the matrix (i.e. ``A.ndim()`` and ``A.shape().size()``), though it is no problem to reference to a matrix certain index using a higher-dimensional equivalent. For example:

.. code-block:: cpp

  cppmat::matrix<double> A({10,10});

  A(5,5,0) = ...

is perfectly acceptable. Note that higher-dimensions can only be trailing ones, using for example ``A(0,5,5)`` is not acceptable, nor is, of course, ``A(5,5,1)``.

Similarly to refer to the beginning of a block (e.g. a row) one can omit the zero arguments. For example, to the beginning of the second row of the above matrix one can use ``&A(1)``.

.. _matrix-index-advanced:

Advanced indexing
-----------------

To allow an arbitrary number of indices at runtime (i.e. the case in which the number of indices is not known at compile time), ``cppmat::matrix`` can also be supplied with the indices stored in a list, using the ``.at(first,last)``. One the just supplied iterators to the beginning and the end of this list. When the indices are also stored in a ``cppmat::matrix`` these iterators can be easily obtained using ``.item(i,j)``. Consider this example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    // example matrix
    // --------------

    cppmat::matrix<size_t> A({2,4});

    A(0,0) =  0; A(0,1) =  1; A(0,2) =  2; A(0,3) =  3;
    A(1,0) = 10; A(1,1) = 11; A(1,2) = 12; A(1,3) = 13;

    // view, based on list of indices
    // ------------------------------

    cppmat::matrix<size_t> index({2,2});

    index(0,0) = 0; index(0,1) = 1;
    index(1,0) = 1; index(1,1) = 2;

    for ( size_t i = 0 ; i < index.shape(0) ; ++i )
      std::cout << A.at(index.item(i), index.item(i)+index.shape(1)) << std::endl;

    return 0;
  }

View
----

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formatting use the more classic C syntax ``A.printf("%16.8e");``

.. _matrix2:

cppmat::matrix2
===============

[:download:`matrix2.h <../src/cppmat/matrix2.h>`, :download:`matrix2.cpp <../src/cppmat/matrix2.cpp>`]

Class for 2-d matrices. For example:

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

[:download:`vector.h <../src/cppmat/vector.h>`, :download:`vector.cpp <../src/cppmat/vector.cpp>`]

Class for 1-d matrices (a.k.a. vectors). For example:

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

  Compared to ``std::vector`` this class is not so much different, with the exception that it provides indexing also with round brackets, and automated printing of entries.

