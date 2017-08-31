
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
