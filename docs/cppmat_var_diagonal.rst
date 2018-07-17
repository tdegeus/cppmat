
****************
cppmat::diagonal
****************

.. _var_diagonal_matrix:

cppmat::diagonal::matrix
========================

[:download:`var_diagonal_matrix.h <../src/cppmat/var_diagonal_matrix.h>`, :download:`var_diagonal_matrix.hpp <../src/cppmat/var_diagonal_matrix.hpp>`]

Square, diagonal, matrix, whereby only the diagonal components are stored:

.. code-block:: none

  [ X       ;
       X    ;
          X ]

*The remaining components are imposed to be zero*. This offers memory advantages, but also computational advantages as the library is fed with additional knowledge of the matrix.

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::diagonal::matrix<double> A(3,3);

      A(0,0) = ...

      // A(0,1) = ... -> not allowed
      // ... = A(0,1) -> allowed, returns zero

      ...

      std::cout << A << std::endl;

      return 0;
  }

Storage
-------

The storage order is as follows:

.. code-block:: none

  [ 0       ;
       1    ;
          2 ]
