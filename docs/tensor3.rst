
.. _tensor3:

******************
<cppmat/tensor3.h>
******************

A specialization of :ref:`tensor`, that takes advantage of the knowledge that the arrays are fixed-size and relatively small. Also several loops are unrolled. All the functionality :ref:`tensor` is available here, whereby the following names are different.

*   Classes:

    -   ``cppmat::tensor4`` -> ``cppmat::tensor3_4``
    -   ``cppmat::tensor2`` -> ``cppmat::tensor3_2``
    -   ``cppmat::tensor2d`` -> ``cppmat::tensor3_2d``
    -   ``cppmat::tensor2s`` -> ``cppmat::tensor3_2s``
    -   ``cppmat::vector`` -> ``cppmat::vector3``

*   Functions:

    - ``cppmat::identity2`` -> ``cppmat::identity3_2``
    - ``cppmat::identity4`` -> ``cppmat::identity3_4``
    - ``cppmat::identity4s`` -> ``cppmat::identity3_4s``
    - ``cppmat::identity4d`` -> ``cppmat::identity3_4d``
    - ``cppmat::identityII`` -> ``cppmat::identity3_II``

Compared to :ref:`tensor` the dimension argument must be omitted everywhere. The same example as before now looks as follows:

.. code-block:: cpp

  #include <cppmat/tensor3.h>

  int main()
  {
      cppmat::tensor3_4<double> A = cppmat::identity3_4();

      ...

      return 0;
  }
