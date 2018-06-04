
************
cppmat::view
************

.. _map_regular_array:

cppmat::view::array
===================

[:download:`map_regular_array.h <../src/cppmat/map_regular_array.h>`, :download:`map_regular_array.hpp <../src/cppmat/map_regular_array.hpp>`]

This class can be used to 'view' a ``const`` external pointer. This can be useful to refer to a part of a bigger array. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> container = cppmat::array<double>::Arange({100,4,2});

      cppmat::view::array<double,2,4,2> view;

      view.setMap(&container(10)); // equivalent to "view.setMap(&container(10,0,0));"

      std::cout << view << std::endl;
  }

This prints:

.. code-block:: bash

  80, 81;
  82, 83;
  84, 85;
  86, 87;

.. warning::

  Since C++ performs garbage collection you should use ``cppmat::view`` with care. You are responsible that pointers do not go out of scope.

.. tip::

  One can also use the ``Map`` constructor instead of the ``setMap`` method:

  .. code-block:: cpp

    using Mat = cppmat::view::matrix<double,4,2>;

    Mat view = Mat::Map(&container(10));

.. note::

  This function cannot make any modification to the view. Its usage is thus somewhat limited. To get a wider functionality use :ref:`fix_regular_array`. For example:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
        cppmat::array<double> container = cppmat::array<double>::Arange({100,4,2});

        cppmat::tiny::array<double,2,4,2> copy;

        copy.setCopy(container.item(10), container.item(10)+copy.size());

        std::cout << copy << std::endl;
    }

  Note that ``copy`` is now a copy. I.e. any modification to ``copy`` will not result in a modification in ``container``.

  Note that the following syntax could also have been used:

  .. code-block:: cpp

    using Mat = cppmat::tiny::matrix<double,4,2>;

    Mat copy = Mat::Copy(container.item(10), container.item(10)+8);

  Or the following:

  .. code-block:: cpp

    using Mat = cppmat::tiny::matrix<double,4,2>;

    Mat copy = Mat::Copy(container.item(10));

  Or the following:

  .. code-block:: cpp

    std::copy(container.item(10), container.item(10)+copy.size(), copy.data());

.. _map_regular_matrix:

cppmat::view::matrix
====================

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::matrix<double,10,10> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_regular_matrix`.

.. _map_regular_vector:

cppmat::view::vector
====================

[:download:`map_regular_vector.h <../src/cppmat/map_regular_vector.h>`, :download:`map_regular_vector.hpp <../src/cppmat/map_regular_vector.hpp>`]

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::vector<double,10> A;

      A.setMap(...)

      ... = A(0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_regular_vector`.

.. _map_symmetric_matrix:

cppmat::view::symmetric::matrix
===============================

[:download:`map_symmetric_matrix.h <../src/cppmat/map_symmetric_matrix.h>`, :download:`map_symmetric_matrix.hpp <../src/cppmat/map_symmetric_matrix.hpp>`]

Class to view a pointer to a fixed size, symmetric, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::symmetric::matrix<double,10,10> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_symmetric_matrix`.

.. _map_diagonal_matrix:

cppmat::view::diagonal::matrix
==============================

[:download:`map_diagonal_matrix.h <../src/cppmat/map_diagonal_matrix.h>`, :download:`map_diagonal_matrix.hpp <../src/cppmat/map_diagonal_matrix.hpp>`]

Class to view a pointer to a fixed size, symmetric, matrices. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::diagonal::matrix<double,10,10> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_diagonal_matrix`.

.. _map_cartesian:

cppmat::view::cartesian
=======================

.. _map_cartesian_tensor4:

``cppmat::view::cartesian::tensor4``
------------------------------------

[:download:`map_cartesian_tensor4.h <../src/cppmat/map_cartesian_tensor4.h>`, :download:`map_cartesian_tensor4.hpp <../src/cppmat/map_cartesian_tensor4.hpp>`]

Class to view a pointer to a fixed size, fourth order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::cartesian::tensor4<double,3> A;

      A.setMap(...)

      ... = A(0,0,0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_cartesian_tensor4`.

.. _map_cartesian_tensor2:

``cppmat::view::cartesian::tensor2``
------------------------------------

[:download:`map_cartesian_tensor2.h <../src/cppmat/map_cartesian_tensor2.h>`, :download:`map_cartesian_tensor2.hpp <../src/cppmat/map_cartesian_tensor2.hpp>`]

Class to view a pointer to a fixed size, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::cartesian::tensor2<double,3> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_cartesian_tensor2`.

.. _map_cartesian_tensor2s:

``cppmat::view::cartesian::tensor2s``
-------------------------------------

[:download:`map_cartesian_tensor2s.h <../src/cppmat/map_cartesian_tensor2s.h>`, :download:`map_cartesian_tensor2s.hpp <../src/cppmat/map_cartesian_tensor2s.hpp>`]

Class to view a pointer to a fixed size, symmetric, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::cartesian::tensor2s<double,3> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_cartesian_tensor2s`.

.. _map_cartesian_tensor2d:

``cppmat::view::cartesian::tensor2d``
-------------------------------------

[:download:`map_cartesian_tensor2d.h <../src/cppmat/map_cartesian_tensor2d.h>`, :download:`map_cartesian_tensor2d.hpp <../src/cppmat/map_cartesian_tensor2d.hpp>`]

Class to view a pointer to a fixed size, diagonal, second order tensors. For a 3-d tensor

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::cartesian::tensor2d<double,3> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_cartesian_tensor2d`.

.. _map_cartesian_vector:

``cppmat::view::cartesian::vector``
-----------------------------------

[:download:`map_cartesian_vector.h <../src/cppmat/map_cartesian_vector.h>`, :download:`map_cartesian_vector.hpp <../src/cppmat/map_cartesian_vector.hpp>`]

Class to view a pointer to a fixed size, vector. For a 3-d vector

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::view::cartesian::vector<double,3> A;

      A.setMap(...)

      ... = A(0,0)

      ...

      return 0;
  }

Most methods are the same as for :ref:`fix_cartesian_vector`.

