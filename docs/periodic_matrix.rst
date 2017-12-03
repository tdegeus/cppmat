
.. _periodic-matrix:

************************
cppmat::periodic::matrix
************************

This class is identical to :ref:`matrix`, with the exception that periodic indices can be used. For example ``-1`` will refer to the last index along that dimension (i.e. ``-1 -> N-1``, with ``N = A.shape(i)``), while ``N`` will refer the the first index (``N -> 0``).

This does require a check and possible modification for each index reference.
