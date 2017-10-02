
.. _tensor2:

******************
<cppmat/tensor2.h>
******************

A specialization of :ref:`tensor`, that takes advantage of the knowledge that the arrays are fixed-size and small. Also several loops are unrolled. All the functionality and names of :ref:`tensor` is transferable to :ref:`tensor2`, whereby the only difference is in the following names:

*   Classes:

    -   ``cppmat::tensor4`` -> ``cppmat::tensor2_4``
    -   ``cppmat::tensor2`` -> ``cppmat::tensor2_2``
    -   ``cppmat::tensor2d`` -> ``cppmat::tensor2_2d``
    -   ``cppmat::tensor2s`` -> ``cppmat::tensor2_2s``
    -   ``cppmat::vector`` -> ``cppmat::vector2``

*   Functions:

    - ``cppmat::identity2`` -> ``cppmat::identity2_2``
    - ``cppmat::identity4`` -> ``cppmat::identity2_4``
    - ``cppmat::identity4s`` -> ``cppmat::identity2_4s``
    - ``cppmat::identity4d`` -> ``cppmat::identity2_4d``
    - ``cppmat::identityII`` -> ``cppmat::identity2_II``

Compared to :ref:`tensor` the dimension argument must be omitted everywhere. The same example as before now looks as follows:

.. code-block:: cpp

  #include <cppmat/tensor2.h>

  int main()
  {
      cppmat::tensor2_4<double> A = cppmat::identity2_4();

      ...

      return 0;
  }

.. note::

  The is no way to automatically switch between :ref:`tensor2` and :ref:`tensor`
