
.. _matrix:

*****************
<cppmat/matrix.h>
*****************

Header-only module that provides a C++ class for n-d matrices. For example a rank 3 matrix is allocated as follows:

.. code-block:: cpp

  #include <cppmat/matrix.h>

  int main()
  {
      cppmat::matrix<double> A({10,10,10});

      A(0,0,0) = ...

      ...

      return 0;
  }

.. note:: **Tip**

  To print, use the common C++ ``std::cout << A << std::endl;``. To customize formating use the more classic C syntax ``A.printf("%16.8e");``

.. note:: **Tip**

  It is no problem to reference to a matrix index using a higher-dimensional index. For example:

  .. code-block:: cpp

    cppmat::matrix<double> A({10,10});

    A(5,5,0) = ...

  is perfectly acceptable. Note that the higher-dimensions should be the trailing ones, e.g. using for example ``A(0,5,5)`` is not acceptable, nor is of course ``A(5,5,1)``.

Methods
=======

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
