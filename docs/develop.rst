
********************
Notes for developers
********************

Make changes / additions
========================

Be sure to run the verification codes in ``develop/``! All existing checks should pass, while new check should be added to check new functionality.

Python
======

The Python package of this module ``cppmat/__init__.py`` is essentially used to allow distribution of the header files that constitute this library through PyPi. In addition a small Python package *cppmat* is provided that allows easy ``setup.py`` formulations of derived packages. These features can also be used when one is just interested in using *pybind11* and one does not intend to use *cppmat* itself.

Create a new release
====================

1.  Update the version numbers as follows:

    -   Modify ``__version__`` in ``setup.py``.
    -   Modify ``version`` in ``cppmat.pc.in``

2.  Upload the changes to GitHub and create a new release there (with the correct version number).

3.  Upload the package to PyPi:

    .. code-block:: bash

      $ python3 setup.py bdist_wheel --universal
      $ twine upload dist/*
