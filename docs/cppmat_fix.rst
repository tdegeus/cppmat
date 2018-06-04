
************
cppmat::tiny
************

.. _fix_regular_array:

cppmat::tiny::array
====================

[:download:`fix_regular_array.h <../src/cppmat/fix_regular_array.h>`, :download:`fix_regular_array.hpp <../src/cppmat/fix_regular_array.hpp>`]

Class for fixed size, small, n-d arrays. For example for a rank 3 array:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::array<double,3,10,10,10> A;

      A(0,0,0) = ...

      ...

      return 0;
  }

Note that the first 'shape' is the rank of the array, the rest are the shape along each axis.

Compared to :ref:`var_regular_array` the size of the array cannot be dynamically changed. Consequently there is no dynamic memory allocation, often resulting in faster behavior. For the rest, most methods are the same as for :ref:`var_regular_array`, though sometimes slightly more limited in use.

.. _fix_regular_matrix:

cppmat::tiny::matrix
====================

[:download:`fix_regular_matrix.h <../src/cppmat/fix_regular_matrix.h>`, :download:`fix_regular_matrix.hpp <../src/cppmat/fix_regular_matrix.hpp>`]

Class for fixed size, small, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::matrix<double,10,10> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_regular_matrix`.

.. _fix_regular_vector:

cppmat::tiny::vector
====================

[:download:`fix_regular_vector.h <../src/cppmat/fix_regular_vector.h>`, :download:`fix_regular_vector.hpp <../src/cppmat/fix_regular_vector.hpp>`]

Class for fixed size, small, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::vector<double,10> A;

      A(0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_regular_vector`.

.. _fix_symmetric_matrix:

cppmat::tiny::symmetric::matrix
===============================

[:download:`fix_symmetric_matrix.h <../src/cppmat/fix_symmetric_matrix.h>`, :download:`fix_symmetric_matrix.hpp <../src/cppmat/fix_symmetric_matrix.hpp>`]

Class for fixed size, small, symmetric, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::symmetric::matrix<double,10,10> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_symmetric_matrix`.

.. _fix_diagonal_matrix:

cppmat::tiny::diagonal::matrix
==============================

[:download:`fix_diagonal_matrix.h <../src/cppmat/fix_diagonal_matrix.h>`, :download:`fix_diagonal_matrix.hpp <../src/cppmat/fix_diagonal_matrix.hpp>`]

Class for fixed size, small, symmetric, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::diagonal::matrix<double,10,10> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_diagonal_matrix`.

.. _fix_cartesian:

cppmat::tiny::cartesian
=======================

.. _fix_cartesian_tensor4:

``cppmat::tiny::cartesian::tensor4``
------------------------------------

[:download:`fix_cartesian_tensor4.h <../src/cppmat/fix_cartesian_tensor4.h>`, :download:`fix_cartesian_tensor4.hpp <../src/cppmat/fix_cartesian_tensor4.hpp>`]

Class for fixed size, small, fourth order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::cartesian::tensor4<double,3> A;

      A(0,0,0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_cartesian_tensor4`.

.. _fix_cartesian_tensor2:

``cppmat::tiny::cartesian::tensor2``
------------------------------------

[:download:`fix_cartesian_tensor2.h <../src/cppmat/fix_cartesian_tensor2.h>`, :download:`fix_cartesian_tensor2.hpp <../src/cppmat/fix_cartesian_tensor2.hpp>`]

Class for fixed size, small, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::cartesian::tensor2<double,3> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_cartesian_tensor2`.

.. _fix_cartesian_tensor2s:

``cppmat::tiny::cartesian::tensor2s``
-------------------------------------

[:download:`fix_cartesian_tensor2s.h <../src/cppmat/fix_cartesian_tensor2s.h>`, :download:`fix_cartesian_tensor2s.hpp <../src/cppmat/fix_cartesian_tensor2s.hpp>`]

Class for fixed size, small, symmetric, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::cartesian::tensor2s<double,3> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_cartesian_tensor2s`.

.. _fix_cartesian_tensor2d:

``cppmat::tiny::cartesian::tensor2d``
-------------------------------------

[:download:`fix_cartesian_tensor2d.h <../src/cppmat/fix_cartesian_tensor2d.h>`, :download:`fix_cartesian_tensor2d.hpp <../src/cppmat/fix_cartesian_tensor2d.hpp>`]

Class for fixed size, small, diagonal, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::cartesian::tensor2d<double,3> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_cartesian_tensor2d`.

.. _fix_cartesian_vector:

``cppmat::tiny::cartesian::vector``
-----------------------------------

[:download:`fix_cartesian_vector.h <../src/cppmat/fix_cartesian_vector.h>`, :download:`fix_cartesian_vector.hpp <../src/cppmat/fix_cartesian_vector.hpp>`]

Class for fixed size, small, vector. For a 3-d vector

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::tiny::cartesian::vector<double,3> A;

      A(0,0) = ...

      ...

      return 0;
  }

Most methods are the same as for :ref:`var_cartesian_vector`.
