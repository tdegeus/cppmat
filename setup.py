desc = '''
Module that provides a header-only library that contains an n-d matrix class and fourth-, second-,
and first-order tensors in C++. Each of these classes has a direct interface to NumPy through
pybind11.
'''

from setuptools import setup

__version__ = '0.1.8'

setup(
   name='cppmat',
   version=__version__,
   description='Matrix (n-d) and tensors in C++',
   author='Tom de Geus',
   author_email='tom@geus.me',
   url='https://github.com/tdegeus/cppmat',
   license='MIT',
   packages=['cppmat'],
   headers=[
     'include/cppmat/matrix.h',
     'include/cppmat/tensor.h',
     'include/cppmat/pybind11_matrix.h',
     'include/cppmat/pybind11_tensor.h',
   ],
   install_requires=['pybind11>=2.1.0'],
   keywords='C++11, Python bindings',
   long_description=desc
)
