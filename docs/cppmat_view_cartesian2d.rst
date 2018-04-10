
.. _view_cartesian2d:

*************************
cppmat::view::cartesian2d
*************************

[:download:`view_tensor2.h <../src/cppmat/view_tensor2.h>`, :download:`view_tensor2.cpp <../src/cppmat/view_tensor2.cpp>`]

This class can be used to 'view' and external pointer. This can be useful to refer to a part of a bigger matrix. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> container({50,50,3});

      container.arange();

      cppmat::view::cartesian2d::tensor2s<double> view;

      view.map(&container(10,0)); // equivalent to "view.map(&container(10,0,0));"

      std::cout << view << std::endl;
  }

This prints:

.. code-block:: bash

  1500, 1501;
  1501, 1502;

.. note::

  This function cannot make any modification to the view. Its usage is thus somewhat limited. To get a wider functionality use :ref:`cartesian2d`. For example:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
        cppmat::matrix<double> container({50,50,3});

        container.arange();

        cppmat::cartesian2d::tensor2s<double> view;

        std::copy(container.item(10,10), container.item(10,10)+view.size(), view.data());

        std::cout << view << std::endl;
    }

  Note that ``view`` is now a copy. I.e. any modification to ``view`` will not result in a modification in ``container``.

.. warning::

  Since C++ performs garbage collection you should use ``cppmat::view`` with care. You are responsible that pointers do not go out of scope.

Classes
=======

.. _view_cartesian2d_tensor4:

cppmat::view::cartesian2d::tensor4
----------------------------------

4th-order tensor (rank 4 tensor).

.. _view_cartesian2d_tensor2:

cppmat::view::cartesian2d::tensor2
----------------------------------

2nd-order tensor (rank 2 tensor).

.. _view_cartesian2d_tensor2s:

cppmat::view::cartesian2d::tensor2s
-----------------------------------

Symmetric 2nd-order tensor.

.. _view_cartesian2d_tensor2d:

cppmat::view::cartesian2d::tensor2d
-----------------------------------

Diagonal 2nd-order tensor.

.. _view_cartesian2d_vector:

cppmat::view::cartesian2d::vector
---------------------------------

Vector (rank 1 tensor).
