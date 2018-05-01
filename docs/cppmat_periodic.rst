
****************
cppmat::periodic
****************

The classes below are identical to those under :ref:`cppmat`, with the exception that periodic indices can be used. For example ``-1`` will refer to the last index along that dimension (i.e. ``-1 -> N-1``, with ``N = A.shape(i)``), while ``N`` will refer to the first index (``N -> 0``).

.. _periodic_array:

cppmat::periodic::array
=======================

[:download:`periodic_array.h <../src/cppmat/periodic_array.h>`, :download:`periodic_array.cpp <../src/cppmat/periodic_array.cpp>`]


.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::array<double> A({10,10,10});

      A(-1,-1,-1) = ... // last item

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`regular_array`.

.. _periodic_matrix:

cppmat::periodic::matrix
========================

[:download:`periodic_matrix.h <../src/cppmat/periodic_matrix.h>`, :download:`periodic_matrix.cpp <../src/cppmat/periodic_matrix.cpp>`]

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A(10,10);

      A(-1,-1) = ... // last item

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`regular_matrix`.

.. _periodic_vector:

cppmat::periodic::vector
========================

[:download:`periodic_vector.h <../src/cppmat/periodic_vector.h>`, :download:`periodic_vector.cpp <../src/cppmat/periodic_vector.cpp>`]

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::vector<double> A(10);

      A(-1) = ... // last item

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`regular_vector`.
