from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

import sys
import setuptools
import pybind11
import cppmat

__version__ = '0.0.1'

ext_modules = [
  Extension(
    'tensorlib',
    ['tensorlib.cpp'],
    include_dirs=[
      pybind11.get_include(False),
      pybind11.get_include(True ),
      cppmat  .get_include(False),
      cppmat  .get_include(True )
    ],
    language='c++'
  ),
]

class BuildExt(build_ext):
  c_opts = {
    'msvc': ['/EHsc'],
    'unix': [],
  }

  if sys.platform == 'darwin':
    c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

  def build_extensions(self):
    ct = self.compiler.compiler_type
    opts = self.c_opts.get(ct, [])
    if ct == 'unix':
      opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
      opts.append(cppmat.cpp_flag(self.compiler))
    elif ct == 'msvc':
      opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
    for ext in self.extensions:
      ext.extra_compile_args = opts
    build_ext.build_extensions(self)

setup(
  name='tensorlib',
  version=__version__,
  author='Tom de Geus',
  author_email='tom@geus.me',
  url='https://github.com/tdegeus/cppmat/docs/examples/tensorlib',
  description='Tensorlib',
  long_description='',
  license='MIT',
  ext_modules=ext_modules,
  install_requires=['pybind11>=2.1.0','cppmat>=0.1.5'],
  cmdclass={'build_ext': BuildExt},
  zip_safe=False,
)
