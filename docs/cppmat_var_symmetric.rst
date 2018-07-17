
*****************
cppmat::symmetric
*****************

.. _var_symmetric_matrix:

cppmat::symmetric::matrix
=========================

[:download:`var_symmetric_matrix.h <../src/cppmat/var_symmetric_matrix.h>`, :download:`var_symmetric_matrix.hpp <../src/cppmat/var_symmetric_matrix.hpp>`]

Square, symmetric, matrix, whereby only the upper-diagonal components are stored:

.. code-block:: none

  [ X, X, X ;
       X, X ;
          X ]

*The remaining components are inferred from symmetry*. This offers memory advantages, but also computational advantages as the library is fed with additional knowledge of the matrix.

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::symmetric::matrix<double> A(3,3);

      A(0,0) = ...

      // A(0,1) = ... -> same as A(1,0) = ...

      ...

      std::cout << A << std::endl;

      return 0;
  }

Storage
-------

The storage order is as follows:

.. code-block:: none

  [ 0, 1, 2 ;
       3, 4 ;
          5 ]

For an ``NxN`` matrix (square by definition), the component ``(i,j)`` can be extracted by

.. code-block:: cpp

  if (i <= j) i*N - (i-1)*i/2 + j - i;
  else        j*N - (j-1)*j/2 + i - j;
