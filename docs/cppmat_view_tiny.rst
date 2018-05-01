
.. _view:

************
cppmat::view
************

.. _view_tiny_matrix:

cppmat::view::matrix
====================

[:download:`view_tiny_matrix.h <../src/cppmat/view_tiny_matrix.h>`, :download:`view_tiny_matrix.cpp <../src/cppmat/view_tiny_matrix.cpp>`]

This class can be used to 'view' and external pointer. This can be useful to refer to a part of a bigger matrix. For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> container = cppmat::array<double>::Arange({100,4,2});

      cppmat::view::matrix<double,4,2> view;

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

.. note::

  One can also use the ``Map`` constructor instead of the ``setMap`` method:

  .. code-block:: cpp

    using Mat = cppmat::view::matrix<double,4,2>;

    Mat view = Mat::Map(&container(10));

.. note::

  This function cannot make any modification to the view. Its usage is thus somewhat limited. To get a wider functionality use :ref:`tiny_matrix`. For example:

  .. code-block:: cpp

    #include <cppmat/cppmat.h>

    int main()
    {
        cppmat::array<double> container = cppmat::array<double>::Arange({100,4,2});

        cppmat::tiny::matrix<double,4,2> copy;

        copy.setCopy(container.item(10), container.item(10)+copy.size());

        std::cout << copy << std::endl;
    }

  Note that ``copy`` is now a copy. I.e. any modification to ``copy`` will not result in a modification in ``container``. Also note that one could have used

  .. code-block:: cpp

    using Mat = cppmat::tiny::matrix<double,4,2>;

    Mat copy = Mat::Copy(container.item(10), container.item(10)+8);

  or

  .. code-block:: cpp

    std::copy(container.item(10), container.item(10)+copy.size(), copy.data());

.. _view_tiny_vector:

cppmat::view::vector
====================

[:download:`view_tiny_vector.h <../src/cppmat/view_tiny_vector.h>`, :download:`view_tiny_vector.cpp <../src/cppmat/view_tiny_vector.cpp>`]

See :ref:`view_tiny_matrix` and :ref:`tiny_vector`.
