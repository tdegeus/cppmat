
.. _tensor2:

******************
cppmat::cartesian2
******************

A specialization of :ref:`tensor`, that takes advantage of the knowledge that the arrays are fixed-size and small. Also several loops are unrolled. All the functionality and names of :ref:`tensor` are transferable to :ref:`tensor2`.

Compared to :ref:`tensor` the dimension argument must be omitted everywhere. The same example as before now looks as follows:

.. code-block:: cpp

  #include <cppmat/tensor2.h>

  int main()
  {
      cppmat::cartesian2d::tensor4<double> A = cppmat::cartesian2d::identity4<double>();

      ...

      return 0;
  }

.. note::

  There is no way to automatically switch between :ref:`tensor2` and :ref:`tensor`.

Map external pointer
--------------------

Like in :ref:`tiny`, the classes under :ref:`tensor2` can be used to 'view' an external pointer. For example, for a matrix of 2-d symmetric tensors:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
      cppmat::matrix2<double> container({50,50,3});

      cppmat::cartesian2d::tensor2s view;

      for ( size_t i = 0 ; i < container.shape(0) ; ++i )
      {
          for ( size_t j = 0 ; j < container.shape(1) ; ++j )
          {
              view.map(&container(i,j));

              view(0,0) = ... // directly stored in "container"
          }
      }
  }

.. _tensor3:

******************
cppmat::cartesian3
******************

Identical to :ref:`tensor2`, but for 3 dimensions.
