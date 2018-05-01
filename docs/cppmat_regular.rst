
.. _cppmat:

******
cppmat
******

.. _regular_array:

cppmat::array
=============

[:download:`regular_array.h <../src/cppmat/regular_array.h>`, :download:`regular_array.cpp <../src/cppmat/regular_array.cpp>`]

Header-only module that provides a C++ class for n-d matrices. For example, a rank 3 matrix is allocated as follows:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

.. note:: **Tip**

  If you know that you will work exclusively with a 2-dimensional matrix, please consider using :ref:`regular_matrix` instead of :ref:`regular_array`. This is generally more efficient. If, on top of that, you want to use only a small matrix, please the fixed sized :ref:`tiny_matrix`. It doesn't need dynamic memory allocation, and can therefore be considerably faster. For the sake of generality there exist also the classes :ref:`regular_vector` and :ref:`tiny_vector`, which have many similarities to the standard C++ ``std::vector``.

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

*   ``A.data()``, ``A.begin()``, ``A.end()``

    Return an iterator to the data, the first, or the last entry of the matrix.

*   ``A.ndim()``

    Returns the number of dimensions of the matrix.

*   ``A.size()``

    Returns the total number of entries in the matrix.

*   ``A.shape(i)``

    Returns the shape along dimension ``i``.

*   ``A.shape()``

    Returns the shape along all dimensions (vector).

*   ``A.resize({...})``

    Resize the matrix.

*   ``A.reshape({...})``

    Change the shape of the matrix. It is required that the total number of entries does not change.

*   ``A.chdim(N)``

    Change the number of dimensions to ``N``. This affects the outputted ``shape``. For example:

    .. code-block:: cpp

      cppmat::array<double> A({10,10});

      A.chdim(3);

    Has the result that ``A.shape() == {10,10,1}``.

*   ``A.setZero()``, ``A.setOnes()``, ``A.setConstant(...)``, ``A.setArange()``

    Set all entries to zero or one, a constant, or the index in the flat storage.

*   ``A.minCoeff()``, ``A.maxCoeff()``

    Return the minimum or the maximum entry.

*   ``A.sum([axis])``

    Return the sum of all entries, or along one or more axes.

*   ``A.mean([axis])``

    Return the mean of all entries, or along one or more axes.

*   ``A.average(weights[, axis][, normalize])``

    Compute the weighted average of all entries, or along one or more axes. See `NumPy <https://docs.scipy.org/doc/numpy/reference/generated/numpy.average.html>`_  and `Wikipedia <https://en.wikipedia.org/wiki/Weighted_arithmetic_mean>`_. Optionally the result can be returned without normalization.

(Named) constructors
--------------------

*   ``cppmat::array<double>(shape)``

    Allocate to a certain shape, nothing is initialized.

*   ``cppmat::array<double>::Arange(shape)``

    Allocate to a certain shape, set entries to its index in the flat storage.

*   ``cppmat::array<double>::Zero(shape)``

    Allocate to a certain shape, set all entries to zero.

*   ``cppmat::array<double>::Ones(shape)``
*
    Allocate to a certain shape, set all entries to one.

*   ``cppmat::array<double>::Constant(shape, constant)``
*
    Allocate to a certain shape, set all entries to a certain constant.

*   ``cppmat::array<double>::Copy(shape, first, last)``
*
    Allocate to a certain shape, copy the individual entries from some external object that is specified using iterators. Note that the flat-size has to match, i.e. ``last - first == size()``.

.. _matrix-index:

Indexing
--------

In principle the number of indices should match the dimensions of the matrix (i.e. ``A.ndim()``). Though one can:

*   Reference to a certain index using a higher-dimensional equivalent. For example:

    .. code-block:: cpp

      cppmat::array<double> A({10,10});

      A(5,5,0) = ...

    is perfectly acceptable. Note that higher-dimensions can only be trailing ones, using for example ``A(0,5,5)`` is not acceptable, nor is, of course, ``A(5,5,1)``.

*   Refer to the beginning of a block (e.g. a row) by omitting the trailing zero indices. For example, a pointer to the beginning of the second row of the above matrix is obtained by ``&A(1)`` (which is fully equivalent to ``&A(1,0)``).

.. _matrix-iterators:

Iterators
---------

One can obtain iterators to:

*   The beginning of the matrix:

    .. code-block:: cpp

      A.begin()

*   The end of the matrix:

    .. code-block:: cpp

      A.end()

*   A specific point in the matrix

    .. code-block:: cpp

      A.item(i,j,k)

*   The data:

    .. code-block:: cpp

      A.data()

View
----

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formatting use the more classic C syntax ``A.printf("%16.8e");``

.. _matrix-index-advanced:

Advanced indexing
-----------------

To allow an arbitrary number of indices at runtime (i.e. the case in which the number of indices is not known at compile time), ``cppmat::array`` can also be supplied with the indices stored in a list, using the ``.at(first,last)``, where ``first`` and ``last`` are iterators to the beginning and the end of this list of indices. When the indices are also stored in a ``cppmat::array`` these iterators can be easily obtained using ``.item(i,j)``. Consider this example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    // example matrix
    // --------------

    cppmat::array<size_t> A({2,4});

    A(0,0) =  0; A(0,1) =  1; A(0,2) =  2; A(0,3) =  3;
    A(1,0) = 10; A(1,1) = 11; A(1,2) = 12; A(1,3) = 13;

    // view, based on list of indices
    // ------------------------------

    cppmat::array<size_t> index({2,2});

    index(0,0) = 0; index(0,1) = 1;
    index(1,0) = 1; index(1,1) = 2;

    for ( size_t i = 0 ; i < index.shape(0) ; ++i )
      std::cout << A.at(index.item(i), index.item(i)+index.shape(1)) << std::endl;

    return 0;
  }

Storage
-------

The matrix is stored `row-major <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_. For a 2-d matrix of size (3,4) this implies the following storage

.. code-block:: python

  [[0, 1, 2, 3],
   [4, 5, 6, 7]]

The ``strides`` indicate per axis how many entries one needs to skip to proceed to the following entry along that axis. For this example

.. code-block:: python

  strides = [4, 1]

.. note:: References

  *   `Row- and column-major order (Wikipedia) <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_
  *   `Reduction (sum) along arbitrary axes of a multidimensional array (StackOverflow) <https://stackoverflow.com/a/49905058/2646505>`_

.. note::

  One can switch back-and-forth between matrix indices and the plain storage using the ``compress`` and ``decompress`` functions. For example:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
      cppmat::array<size_t> A({2,4});

      std::cout << A.compress(1,2) << std::endl;

      std::vector<size_t> idx = A.decompress(6);

      for ( auto &i : idx )
        std::cout << i << ", ";
      std::cout << std::endl;

      return 0;
    }

  Prints

  .. code-block:: python

    6
    1, 2,

.. _regular_matrix:

cppmat::matrix
==============

[:download:`regular_matrix.h <../src/cppmat/regular_matrix.h>`, :download:`regular_matrix.cpp <../src/cppmat/regular_matrix.cpp>`]

Class for 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A(10,10);

      A(0,0) = ...

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`regular_array`, though there is obviously no ``chdim`` method.

.. _regular_vector:

cppmat::vector
==============

[:download:`regular_vector.h <../src/cppmat/regular_vector.h>`, :download:`regular_vector.cpp <../src/cppmat/regular_vector.cpp>`]

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

  The entire interface is the same as for :ref:`regular_array`, though there is obviously no ``chdim`` method.

.. note::

  Compared to ``std::vector`` this class is not much difference. The only exception that it provides indexing also with round brackets, and automated printing of entries.

