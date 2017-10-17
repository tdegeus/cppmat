
.. _vector:

******************
<cppmat/vector.h>
******************

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

View
====

To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``
