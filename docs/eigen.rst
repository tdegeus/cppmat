
******************
Interface to Eigen
******************

Where possible an interface to the popular Eigen linear algebra library is provided. This practically means that you can copy seamlessly to and from Eigen matrices, but you can also map an Eigen array to one of cppmat's classes (provided that the row-major) storage is used. To enable this feature:

*   Include Eigen before cppmat:

    .. code-block:: cpp

      #include <Eigen/Eigen>
      #include <cppmat/cppmat.h>

*   Define ``CPPMAT_EIGEN`` somewhere before including cppmat:


    .. code-block:: cpp

      #define CPPMAT_EIGEN
      #include <cppmat/cppmat.h>
      #include <Eigen/Eigen>

For example:

.. code-block:: cpp

  #include <Eigen/Eigen>
  #include <cppmat/cppmat.h>

  int main()
  {

    cppmat::cartesian2d::tensor2s<double> A;

    A(0,0) = 10.;
    A(0,1) =  1.;
    A(1,1) = 20.;

    Eigen::Matrix<double,2,2,Eigen::RowMajor> B = A;

    std::cout << "A = " << std::endl << A << std::endl;
    std::cout << "B = " << std::endl << B << std::endl;

  }

will result in

.. code-block:: bash

  A =
  10, 1;
  1, 20;

  B =
  10  1
   1 20
