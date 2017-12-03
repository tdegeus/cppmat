
*****************
Access to storage
*****************

The storage of all the classes can be accessed through the ``data()`` method, which is complemented with the iterators ``begin()`` and ``end()``. Consider the following example:

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  using T2  = cppmat::cartesian::tensor2 <double>;
  using T2s = cppmat::cartesian::tensor2s<double>;
  using T2d = cppmat::cartesian::tensor2d<double>;
  using V   = cppmat::cartesian::vector  <double>;

  int main()
  {

    T2 A(3);

    std::vector<double> data(3*3);

    for ( size_t i = 0 ; i < 3*3 ; ++i )
      data[i] = static_cast<double>(i);

    std::copy(data.begin(), data.end(), A.begin());

    std::cout << "A = " << A << std::endl;

    T2s B = A.cast<T2s>();

    std::cout << "B = " << B << std::endl;

    V C(3);

    C[0] = 10.01;
    C[1] = 20.01;
    C[2] = 30.01;

    std::cout << "C = " << std::setw(10) << std::setprecision(6) << C << std::endl;

    std::copy(C.begin(), C.end(), &A(1,0));

    std::cout << "converted A = " << A << std::endl;

  }

