
******
cppmat
******

.. _var_regular_array:

cppmat::array
=============

[:download:`var_regular_array.h <../src/cppmat/var_regular_array.h>`, :download:`var_regular_array.hpp <../src/cppmat/var_regular_array.hpp>`]

A C++ class for dynamically sized arrays or arbitrary rank. For example, a rank 3 array is allocated as follows:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      std::cout << A << std::endl;

      return 0;
  }

.. tip::

  *  If you know that you will work exclusively with a rank 1 or 2 array (i.e. a vector or a matrix), consider using :ref:`var_regular_vector`, :ref:`var_regular_matrix`, :ref:`var_symmetric_matrix`, and , :ref:`var_diagonal_matrix`. This can enhance readability and/or efficiency.

  *  If your array is not very big and its size is known at compile time consider using :ref:`fix_regular_array` (or the fixed size equivalents of the other classes). This avoids dynamic memory allocation, and usually speeds-up your code.

  *  If your array is part of an external array (for example a bigger array) which you want to just read from, consider using :ref:`map_regular_array`.

  *  To format the print use the regular C++ mechanism, e.g. ``std::cout << std::setw(5) << std::setprecision(3) << A << std::endl;``

Methods
-------

*   ``A(i,j,k)``

    Returns the entry at ``(i,j,k)``. Use this to read or write.

    A negative index may also be used (in that case the indices have to be ``int``) which counts down from the last index along that axis. For example ``A(-1,-1,-1)`` in the last index of the above array. This implies some extra operations, so if you do not use this feature input the indices as ``size_t``.

    The number of indices (i.e. ``A(i)``, ``A(i,j)``, ``A(i,j,k)``, ...) may be lower or equal to the rank, all 'omitted' indices are assumed to be zero.

    See :ref:`array-index` for additional directives.

*   ``A[i]``

    Returns the ``i``-th entry of the plain storage. Use this to read or write.

*   ``A.at(first, last)``

    Returns the entry ``{i,j,k}``, which are stored in a list. The function takes an iterator to the first and the last index of this list. See :ref:`array-index-advanced`.

*   ``A.item(i,j,k)``

    Returns an iterator to the entry at ``(i,j,k)``.

*   ``A.index(i)``

    Returns an iterator to the ``i``-th entry of the plain storage.

*   ``A.data()``, ``A.begin()``, ``A.end()``

    Return an iterator to the data, the first, or the last entry of the matrix.

*   ``A.rank()``

    Returns the ranks of the array (i.e. the number of axes).

*   ``A.size()``

    Returns the total number of entries in the matrix.

*   ``A.shape(i)``

    Returns the shape along dimension ``i`` (a negative number may be used that counts down from the last axis, e.g. ``A.shape(-1)`` is the same as ``A.shape(A.rank()-1)``.

*   ``A.shape()``

    Returns the shape along all dimensions (vector).

*   ``A.resize({...}[, D])``

    Resize the matrix. Enter a value to initialize all allocated entries.

*   ``A.reshape({...})``

    Change the shape of the matrix. It is required that the total number of entries does not change.

*   ``A.chrank(N)``

    Change the rank to ``N`` (with shape 1 along the added axes). A reduction of rank is only allowed if the shape is 1 along the reduced axes.

*   ``A.setZero()``, ``A.setOnes()``, ``A.setConstant(D)``, ``A.setArange()``, ``A.setRandom([start, end])``

    Set all entries to zero or one, a constant, the index in the flat storage, or a random value.

*   ``A.setCopy(first[, last])``

    Copy the individual entries from some external object that is specified using iterators. Note that the flat-size has to match, i.e. ``last - first == size()``.

*   ``A.copyTo(first[, last])``

    Copy the individual entries to an external iterator.

*   ``A.abs()``

    Returns an array with the absolute values of each entry.

*   ``A.norm()``

    Returns the norm (sum of absolute values).

*   ``A.argmin()``, ``A.argmax()``

    Return the plain storage index of the minimum/maximum.

*   ``A.min([axis])``, ``A.max([axis])``

    Return the minimum or the maximum entry.

*   ``A.sum([axis])``

    Return the sum of all entries, or along one or more axes.

*   ``A.mean([axis])``

    Return the mean of all entries, or along one or more axes.

*   ``A.average(weights[, axis, normalize])``

    Compute the weighted average of all entries, or along one or more axes. See `NumPy <https://docs.scipy.org/doc/numpy/reference/generated/numpy.average.html>`_  and `Wikipedia <https://en.wikipedia.org/wiki/Weighted_arithmetic_mean>`_. Optionally the result can be returned without normalization.

*   ``A.where()``

    Returns a vector with the plain storage indices of all non-zero entries.

*   ``A.equal(D)``, ``A.not_equal(D)``, ``A.greater(D)``, ``A.greater_equal(D)``, ``A.less(D)``, ``A.less_equal(D)``

    Return array of booleans, based on the condition.

*    ``A.slice(...)``

     Returns a slice of the array. The input are ``std::vector<size_t>`` with the indices to select along that axis (these vectors can be also input using the ``{...}`` syntax). An empty vector (or simply ``{}``) implies that all indices along that axis are selected.

.. tip::

  If you use something other than ``size_t`` as the type for indices (e.g. ``int``), the functions ``size``, ``shape``, ``rank``, and ``strides`` can be templated to directly get the type you want. For example:

  .. code-block:: cpp

    cppmat::array<double> A({10,10,10});

    for ( int i = 0 ; i < A.size<int>() ; ++i )
       ...

(Named) constructors
--------------------

*   ``cppmat::array<double>(shape)``

    Allocate to a certain shape, nothing is initialized. The ``shape`` has to be specified as a ``std::vector<size_t>``, from which the rank is automatically deduced. Alternatively the ``{...}`` notation can be used, to avoid a separate variable.

*   ``cppmat::array<double>::Random(shape[, start, end])``

    Allocate to a certain shape, set entries to a random value.

*   ``cppmat::array<double>::Arange(shape)``

    Allocate to a certain shape, set entries to its index in the flat storage.

*   ``cppmat::array<double>::Zero(shape)``

    Allocate to a certain shape, set all entries to zero.

*   ``cppmat::array<double>::Ones(shape)``

    Allocate to a certain shape, set all entries to one.

*   ``cppmat::array<double>::Constant(shape, constant)``

    Allocate to a certain shape, set all entries to a certain constant.

*   ``cppmat::array<double>::Copy(shape, first[, last])``

    Allocate to a certain shape, copy the individual entries from some external object that is specified using iterators. Note that the flat-size has to match, i.e. ``last - first == size()``.

External operations
-------------------

*   ``cppmat::array<double> = cppmat::min(A, B)``

    Construct an array taking the minimum of two arrays for each entry.

*   ``cppmat::array<double> = cppmat::max(A, B)``

    Construct an array taking the maximum of two arrays for each entry.

.. _array-index:

Indexing
--------

In principle the number of indices should match the rank of the array (i.e. ``A.rank()``). Though one can:

*   Reference to a certain index using a higher-dimensional equivalent. For example:

    .. code-block:: cpp

      cppmat::array<double> A({10,10});

      A(5,5,0) = ...

    is perfectly acceptable. Note that higher-dimensions can only be trailing ones, using for example ``A(0,5,5)`` is not acceptable, nor is, of course, ``A(5,5,1)``.

*   Refer to the beginning of a block (e.g. a row) by omitting the trailing zero indices. For example, a pointer to the beginning of the second row of the above matrix is obtained by ``&A(1)`` (which is fully equivalent to ``&A(1,0)``).

.. tip::

  A negative index may also be used (in that case the indices have to be ``int``) which counts down from the last index along that axis. For example ``A(-1,-1)`` in the last index of the above matrix. To input any *periodic* index (i.e. to turn-off the bound-checks) use ``.setPeriodic(true)`` on the array object. In that case ``A(-1,-1) == A(10,10)`` for the above matrix.

  This does involve some extra operations, so if you do not use this feature input the indices as ``size_t``.

.. _array-index-advanced:

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

.. tip::

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

.. _var_regular_matrix:

cppmat::matrix
==============

[:download:`var_regular_matrix.h <../src/cppmat/var_regular_matrix.h>`, :download:`var_regular_matrix.hpp <../src/cppmat/var_regular_matrix.hpp>`]

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

The entire interface is the same as for :ref:`var_regular_array`, though there is obviously no ``chrank`` method.

.. _var_regular_vector:

cppmat::vector
==============

[:download:`var_regular_vector.h <../src/cppmat/var_regular_vector.h>`, :download:`var_regular_vector.hpp <../src/cppmat/var_regular_vector.hpp>`]

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

The entire interface is the same as for :ref:`var_regular_array`, though there is obviously no ``chrank`` method.

.. tip::

  One can almost seamlessly switch between ``std::vector`` and ``cppmat::vector``. For example the following would work:

  .. code-block:: cpp

    std::vector<double> A = cppmat::vector<double>::Random(10);
