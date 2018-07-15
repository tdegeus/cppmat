
*********************
STL(-like) extensions
*********************

Out-of-place functions
======================

*   **min/max**

    ``cppmat::array<double> C = cppmat::min(A, B);``

std::vector
===========

*   **Formatted print**

    If cppmat is loaded one can also view STL-vectors as easily as

    ``std::cout << A << std::endl;``

*   **Delete item**

    ``A = cppmat::del(A, -1);``

*   **argsort**

    ``std::vector<size_t> idx = cppmat::argsort(A);``

*   **linspace**

    ``std::vector<double> A = cppmat::linspace(0.0, 1.0, 11);``

*   **min/max**

    ``std::vector<double> C = cppmat::min(A, B);``



