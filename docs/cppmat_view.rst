
.. _view:

************
cppmat::view
************

.. _view_matrix2:

cppmat::view::matrix2
=====================

[:download:`view_matrix2.h <../src/cppmat/view_matrix2.h>`, :download:`view_matrix2.cpp <../src/cppmat/view_matrix2.cpp>`]

This class can be used to 'view' and external pointer. This can be useful to refer to a part of a bigger matrix. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> container({100,4,2});

      container.arange();

      cppmat::view::matrix2<double,4,2> view;

      view.map(&container(10)); // equivalent to "view.map(&container(10,0,0));"

      std::cout << view << std::endl;
  }

This prints:

.. code-block:: bash

  80, 81;
  82, 83;
  84, 85;
  86, 87;

.. note::

  This function cannot make any modification to the view. Its usage is thus somewhat limited. To get a wider functionality use :ref:`tiny_matrix2`. For example:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
        cppmat::matrix<double> container({100,4,2});

        container.arange();

        cppmat::tiny::matrix2<double,4,2> view;

        std::copy(container.item(10), container.item(10)+view.size(), view.data());

        std::cout << view << std::endl;
    }

  Note that ``view`` is now a copy. I.e. any modification to ``view`` will not result in a modification in ``container``.

.. warning::

  Since C++ performs garbage collection you should use ``cppmat::view`` with care. You are responsible that pointers do not go out of scope.

.. _view_vector:

cppmat::view::vector
====================

[:download:`view_vector.h <../src/cppmat/view_vector.h>`, :download:`view_vector.cpp <../src/cppmat/view_vector.cpp>`]

See :ref:`view_matrix2`.
