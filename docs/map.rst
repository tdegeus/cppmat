
*********************************
Map existing data without copying
*********************************

It is possible to map the small, fixed size, matrices and tensors to existing data. Therefore each of the following classes is equipped with ``map(X* ptr)`` function:

*   ``cppmat::tiny::matrix2<T,M,N>``
*   ``cppmat::tiny::vector<T,N>``
*   ``cppmat::cartesian2d::tensor4<T>``
*   ``cppmat::cartesian2d::tensor2<T>``
*   ``cppmat::cartesian2d::tensor2s<T>``
*   ``cppmat::cartesian2d::tensor2d<T>``
*   ``cppmat::cartesian2d::vector<T>``
*   ``cppmat::cartesian3d::tensor4<T>``
*   ``cppmat::cartesian3d::tensor2<T>``
*   ``cppmat::cartesian3d::tensor2s<T>``
*   ``cppmat::cartesian3d::tensor2d<T>``
*   ``cppmat::cartesian3d::vector<T>``

All of the usage and interaction stays the same. Of course, the ``=`` operator will still return a copy.

For example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {

    double data[3];

    data[0] =  0.;
    data[1] =  1.;
    data[2] = 20.;

    cppmat::cartesian2d::tensor2s<double> A;

    A.map(data);

    cppmat::cartesian2d::tensor2s<double> B = A;

    // modifies "A", because it points to "data", but not "B" as it is a copy
    data[0] = 10.;

    // modifies only "B"
    B[0] = -1.;

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;

  }

will result in

.. code-block:: bash

  A =
  10, 1;
  1, 20;

  B =
  -1, 1;
  1, 20;
