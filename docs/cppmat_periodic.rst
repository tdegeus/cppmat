
****************
cppmat::periodic
****************

The class below are identical to those under :ref:`cppmat`, with the exception that periodic indices can be used. For example ``-1`` will refer to the last index along that dimension (i.e. ``-1 -> N-1``, with ``N = A.shape(i)``), while ``N`` will refer the the first index (``N -> 0``).

.. _periodic_matrix:

cppmat::periodic::matrix
========================

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix<double> A({10,10,10});

      A(-1,-1,-1) = ... // last item

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`matrix`.

.. _periodic_matrix2:

cppmat::periodic::matrix2
=========================

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix2<double> A(10,10);

      A(-1,-1) = ... // last item

      ...

      return 0;
  }

.. note::

  The entire interface is the same as for :ref:`matrix2`.

.. _periodic_vector:

cppmat::periodic::vector
========================

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

  The entire interface is the same as for :ref:`vector`.
