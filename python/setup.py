desc = '''
Module that provides a header-only library that contains an n-d matrix class and fourth-, second-,
and first-order tensors (a.k.a. vectors) in C++. Each of these classes has a direct interface to
NumPy through pybind11.

The purpose of this Python package is to simplify the installation and distribution of modules that
depend on cppmat. There is no direct use for Python.
'''

import re
from setuptools import setup

header = open('../src/cppmat/macros.h','r').read()
world  = re.split('(.*)(\#define CPPMAT_WORLD_VERSION\ )([0-9]+)(.*)',header)[3]
major  = re.split('(.*)(\#define CPPMAT_MAJOR_VERSION\ )([0-9]+)(.*)',header)[3]
minor  = re.split('(.*)(\#define CPPMAT_MINOR_VERSION\ )([0-9]+)(.*)',header)[3]

__version__ = '.'.join([world,major,minor])

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
      '../src/cppmat/cppmat.h',
      '../src/cppmat/macros.h',
      '../src/cppmat/matrix.h',
      '../src/cppmat/matrix2.h',
      '../src/cppmat/pybind11_cppmat.h',
      '../src/cppmat/pybind11_matrix.h',
      '../src/cppmat/pybind11_matrix2.h',
      '../src/cppmat/pybind11_tensor.h',
      '../src/cppmat/pybind11_tensor2.h',
      '../src/cppmat/pybind11_tensor3.h',
      '../src/cppmat/pybind11_tiny_matrix2.h',
      '../src/cppmat/pybind11_tiny_vector.h',
      '../src/cppmat/pybind11_vector.h',
      '../src/cppmat/tensor.h',
      '../src/cppmat/tensor2.h',
      '../src/cppmat/tensor3.h',
      '../src/cppmat/tiny_matrix2.h',
      '../src/cppmat/tiny_vector.h',
      '../src/cppmat/vector.h'
   ],
   install_requires = ['pybind11>=2.2.0'],
)
