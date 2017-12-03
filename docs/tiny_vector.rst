
.. _tiny_vector:

********************
cppmat::tiny::vector
********************

Header-only module that provides a C++ class for fixed size, small, 1-d matrices (a.k.a. vectors). For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::vector<double,10> A;

      A(0) = ...

      ...

      return 0;
  }

.. note::

  Compared to `std::vector` this class avoid dynamic memory allocation, and is generally faster. Also it provides indexing also with round brackets, and and automated printing of entries.

Methods
=======

*   ``A(i)``

    Return the entry at ``(i)``.

*   ``A[i]``

    Returns the ``i``-th entry in the matrix viewed as a list.

*   ``size_t n = A.shape(i)``

    Returns the shape along dimension ``i`` (where ``i`` can only be ``0`` in this case).

*   ``std::vector<size_t> shape = A.shape()``

    Returns the shape as an array of length one.

.. note::

  This class does not allow for periodic indices, the never loose efficiency.

View
====

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``
