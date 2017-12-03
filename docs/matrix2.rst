
.. _matrix2:

***************
cppmat::matrix2
***************

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

Methods
=======

*   ``A(i,j)``

    Return the entry at ``(i,j)``.

*   ``A[i]``

    Returns the ``i``-th entry in the matrix viewed as a list.

*   ``size_t n = A.shape(i)``

    Returns the shape along dimension ``i``.

*   ``std::vector<size_t> shape = A.shape()``

    Returns the shape along both dimensions.

*   ``A.reshape(m,n)``

    Change the shape of the matrix. It is required that the total number of entries does not change.

View
====

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``
