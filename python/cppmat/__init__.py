
from setuptools.command.build_ext import build_ext
import sys

# --------------------------------------------------------------------------------------------------

def has_flag(compiler, flagname):
  import tempfile
  with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
    f.write('int main (int argc, char **argv) { return 0; }')
    try:
      compiler.compile([f.name], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
      return False
  return True

# --------------------------------------------------------------------------------------------------

def cpp_flag(compiler):
  if   has_flag(compiler,'-std=c++14'): return '-std=c++14'
  elif has_flag(compiler,'-std=c++11'): return '-std=c++11'
  raise RuntimeError('Unsupported compiler: at least C++11 support is needed')

# --------------------------------------------------------------------------------------------------

def get_include(*args,**kwargs):
  import os
  try:
    from pip import locations
    return os.path.dirname(locations.distutils_scheme('cppmat',*args,**kwargs)['headers'])
  except ImportError:
    return 'include'

# --------------------------------------------------------------------------------------------------

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
      opts.append(cpp_flag(self.compiler))
    elif ct == 'msvc':
      opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
    for ext in self.extensions:
      ext.extra_compile_args = opts
    build_ext.build_extensions(self)

# --------------------------------------------------------------------------------------------------

def find_eigen(hint=None):

  # search with pkgconfig
  # ---------------------

  try:

    import pkgconfig

    if pkgconfig.installed('eigen3','>3.0.0'):
      return pkgconfig.parse('eigen3')['include_dirs'][0]

  except:
    pass

  # manual search
  # -------------

  import os,re

  search_dirs = [] if hint is None else hint
  search_dirs += [
    "/usr/local/include/eigen3",
    "/usr/local/homebrew/include/eigen3",
    "/opt/local/var/macports/software/eigen3",
    "/opt/local/include/eigen3",
    "/usr/include/eigen3",
    "/usr/include/local",
    "/usr/include",
  ]

  for d in search_dirs:
    path = os.path.join(d, "Eigen", "Dense")
    if os.path.exists(path):
      vf = os.path.join(d, "Eigen", "src", "Core", "util", "Macros.h")
      if not os.path.exists(vf):
        continue
      src = open(vf, "r").read()
      v1 = re.findall("#define EIGEN_WORLD_VERSION (.+)", src)
      v2 = re.findall("#define EIGEN_MAJOR_VERSION (.+)", src)
      v3 = re.findall("#define EIGEN_MINOR_VERSION (.+)", src)
      if not len(v1) or not len(v2) or not len(v3):
        continue
      v = "{0}.{1}.{2}".format(v1[0], v2[0], v3[0])
      print("Found Eigen version {0} in: {1}".format(v, d))
      return d

  return None
