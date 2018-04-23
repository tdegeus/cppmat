
*****************
Access to storage
*****************

The storage of all the classes can be accessed through the ``data()`` method, which is complemented with the iterators ``begin()`` and ``end()``. Consider the following examples

Creating a cppmat-object
========================

Copy constructor
----------------

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  using T2 = cppmat::cartesian3d::tensor2<double>;

  int main()
  {
    std::vector<double> data(3*3);

    for ( size_t i = 0 ; i < 3*3 ; ++i )
      data[i] = static_cast<double>(i);

    T2 A = T2::Copy(data.begin(), data.end());

    std::cout << "A = " << A << std::endl;
  }

std::copy
---------

.. code-block:: cpp

  #include <cppmat/cppmat.h>

  using T2 = cppmat::cartesian3d::tensor2<double>;

  int main()
  {
    std::vector<double> data(3*3);

    for ( size_t i = 0 ; i < 3*3 ; ++i )
      data[i] = static_cast<double>(i);

    T2 A;

    std::copy(data.begin(), data.end(), A.begin());

    std::cout << "A = " << A << std::endl;
  }
