
.. _tiny_matrix2:

***********************
<cppmat/tiny_matrix2.h>
***********************

Header-only module that provides a C++ class for fixed size, small, 2-d matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::matrix2<double,10,10> A;

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

.. note::

  This class does not allow for periodic indices, the never loose efficiency.

View
====

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``
