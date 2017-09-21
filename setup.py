desc = '''
Module that provides a header-only library that contains an n-d matrix class and fourth-, second-,
and first-order tensors (a.k.a. vectors) in C++. Each of these classes has a direct interface to
NumPy through pybind11.

The purpose of this Python package is to simplify the installation and distribution of modules that
depend on cppmat. There is no direct use for Python.
'''

from setuptools import setup

__version__ = '0.2.9'

setup(
   name             = 'cppmat',
   description      = 'n-D matrix and tensors in C++',
   long_description = desc,
   keywords         = 'C++, C++11, Python bindings, pybind11',
   version          = __version__,
   license          = 'MIT',
   author           = 'Tom de Geus',
   author_email     = 'tom@geus.me',
   url              = 'https://github.com/tdegeus/cppmat',
   packages         = ['cppmat'],
   headers          = [
     'include/cppmat/matrix.h',
     'include/cppmat/tensor.h',
     'include/cppmat/tensor3.h',
     'include/cppmat/pybind11_matrix.h',
     'include/cppmat/pybind11_tensor3.h',
   ],
   install_requires = ['pybind11>=2.1.0'],
)
