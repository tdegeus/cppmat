
.. _tensor3:

********************
cppmat::cartesian3::
********************

A specialization of :ref:`tensor`, that takes advantage of the knowledge that the arrays are fixed-size and small. Also several loops are unrolled. All the functionality and names of :ref:`tensor` are transferable to :ref:`tensor3`.

Compared to :ref:`tensor` the dimension argument must be omitted everywhere. The same example as before now looks as follows:

.. code-block:: cpp

  #include <cppmat/tensor3.h>

  int main()
  {
      cppmat::cartesian3d::tensor4<double> A = cppmat::cartesian3d::identity4();

      ...

      return 0;
  }

.. note::

  There is no way to automatically switch between :ref:`tensor3` and :ref:`tensor`.
