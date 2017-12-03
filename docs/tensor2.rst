
.. _tensor2:

********************
cppmat::cartesian2::
********************

A specialization of :ref:`tensor`, that takes advantage of the knowledge that the arrays are fixed-size and small. Also several loops are unrolled. All the functionality and names of :ref:`tensor` are transferable to :ref:`tensor2`.

Compared to :ref:`tensor` the dimension argument must be omitted everywhere. The same example as before now looks as follows:

.. code-block:: cpp

  #include <cppmat/tensor2.h>

  int main()
  {
      cppmat::cartesian2d::tensor4<double> A = cppmat::cartesian2d::identity4();

      ...

      return 0;
  }

.. note::

  There is no way to automatically switch between :ref:`tensor2` and :ref:`tensor`.
