
*********
Histogram
*********

histogram
---------

.. code-block:: cpp

  template<typename X>
  std::tuple<std::vector<double>, std::vector<double>> histogram(
    const std::vector<X> &data, size_t bins=10, bool density=false, bool return_edges=false
  )

Create a histogram. Returns ``std::tie(P, x)``: the count and the locations on the bins (their midpoints, or their edges if ``return_edges=true``).

histogram_uniform
-----------------

.. code-block:: cpp

  template<typename X>
  std::tuple<std::vector<double>, std::vector<double>> histogram_uniform(
    const std::vector<X> &data, size_t bins=10, bool density=false, bool return_edges=false
  )

Create a histogram such that each bins contains the same number of entries. Returns ``std::tie(P, x)``: the count and the locations on the bins (their midpoints, or their edges if ``return_edges=true``).
