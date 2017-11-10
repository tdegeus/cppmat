from setuptools import setup, Extension

import sys
import setuptools
import pybind11
import cppmat

__version__ = '0.0.1'

ext_modules = [
  Extension(
    'matrixlib',
    ['matrixlib.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True )
    ],
    language='c++'
  ),
]

setup(
  name               = 'matrixlib',
  description        = 'Matrix library',
  long_description   = 'This is an example module, it no real use!',
  keywords           = 'Example, C++, C++11, Python bindings, pybind11',
  version            = __version__,
  license            = 'MIT',
  author             = 'Tom de Geus',
  author_email       = 'tom@geus.me',
  url                = 'https://github.com/tdegeus/cppmat/docs/examples/matrixlib',
  ext_modules        = ext_modules,
  extra_compile_args = ["-DNDEBUG"], # switch off assertions
  install_requires   = ['pybind11>=2.1.0','cppmat>=0.3.0'],
  cmdclass           = {'build_ext': cppmat.BuildExt},
  zip_safe           = False,
)
