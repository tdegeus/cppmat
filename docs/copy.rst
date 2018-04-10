
*****************
Access to storage
*****************

The storage of all the classes can be accessed through the ``data()`` method, which is complemented with the iterators ``begin()`` and ``end()``. Consider the following examples

Example 1: std::copy
====================

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    cppmat::cartesian::tensor2<double> A(3);

    std::vector<double> data(3*3);

    for ( size_t i = 0 ; i < 3*3 ; ++i )
      data[i] = static_cast<double>(i);

    std::copy(data.begin(), data.end(), A.begin());

    std::cout << "A = " << A << std::endl;
  }

Example 2: construct + copy from iterator
=========================================

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  int main()
  {
    std::vector<double> data(3*3);

    for ( size_t i = 0 ; i < 3*3 ; ++i )
      data[i] = static_cast<double>(i);

    cppmat::cartesian::tensor2<double> A(3, data.begin(), data.end());

    std::cout << "A = " << A << std::endl;
  }
