
import os, sys, re

from setuptools                   import Extension
from setuptools.command.build_ext import build_ext

# ==================================================================================================

def has_flag(compiler, flagname):
  r'''
Check if a compiler supports a certain flag. Returns ``True`` or ``False``.

The function creates a temporary file and tries compiling with the compiler flag.

(c) Sylvain Corlay, https://github.com/pybind/python_example
  '''

  import tempfile

  with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:

    f.write('int main (int argc, char **argv) { return 0; }')

    try:
      compiler.compile([f.name], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
      return False

  return True

# ==================================================================================================

def cpp_flag(compiler):
  r'''
Set the c++14 standard, or else the c++11 standard. In case neither can be set the function
raises an error.

(c) Sylvain Corlay, https://github.com/pybind/python_example
  '''

  if   has_flag(compiler,'-std=c++14'): return '-std=c++14'
  elif has_flag(compiler,'-std=c++11'): return '-std=c++11'
  raise RuntimeError('Unsupported compiler: at least C++11 support is needed')

# ==================================================================================================

def get_include(user=False):
  r'''
Get the relevant ``include`` directory.

(c) Sylvain Corlay, https://github.com/pybind/python_example
  '''

  from distutils.dist import Distribution

  # Are we running in a virtual environment?
  # - check
  virtualenv = hasattr(sys, 'real_prefix') or sys.prefix != getattr(sys, "base_prefix", sys.prefix)
  # - return path
  if virtualenv:
    return os.path.join(sys.prefix, 'include', 'site', 'python' + sys.version[:3])

  # Search
  dist = Distribution({'name': 'cppmat'})
  dist.parse_config_files()
  dist_cobj = dist.get_command_obj('install', create=True)

  # Search for packages in user's home directory?
  if user:
    dist_cobj.user = user
    dist_cobj.prefix = ""

  # Search
  dist_cobj.finalize_options()

  return os.path.dirname(dist_cobj.install_headers)

# ==================================================================================================

class BuildExt(build_ext):
  r'''
Define class to build the extension.

(c) Sylvain Corlay, https://github.com/pybind/python_example
  '''

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

# ==================================================================================================

def find_eigen(hint=None):
  r'''
Try to find the Eigen library. If successful the include directory is returned.
  '''

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

# ==================================================================================================

class CMakeExtension(Extension):
  r'''
(c) Dean Moldovan, https://github.com/pybind/cmake_example
  '''

  def __init__(self, name, sourcedir=''):

    Extension.__init__(self, name, sources=[])

    self.sourcedir = os.path.abspath(sourcedir)

# ==================================================================================================

class CMakeBuild(build_ext):
  r'''
(c) Dean Moldovan, https://github.com/pybind/cmake_example
  '''

  # ------------------------------------------------------------------------------------------------

  def run(self):

    import platform, subprocess

    from distutils.version import LooseVersion

    try:
      out = subprocess.check_output(['cmake', '--version'])
    except OSError:
      raise RuntimeError("CMake must be installed to build the following extensions: " +
                           ", ".join(e.name for e in self.extensions))

    if platform.system() == "Windows":

      cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))

      if cmake_version < '3.1.0':
        raise RuntimeError("CMake >= 3.1.0 is required on Windows")

    for ext in self.extensions:

      self.build_extension(ext)

  # ------------------------------------------------------------------------------------------------

  def build_extension(self, ext):

    import platform, subprocess

    extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

    cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                  '-DPYTHON_EXECUTABLE=' + sys.executable]

    cfg = 'Debug' if self.debug else 'Release'

    build_args = ['--config', cfg]

    if platform.system() == "Windows":
      cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
      if sys.maxsize > 2**32:
        cmake_args += ['-A', 'x64']
      build_args += ['--', '/m']
    else:
      cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
      build_args += ['--', '-j2']

    env = os.environ.copy()
    env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                          self.distribution.get_version())

    if not os.path.exists(self.build_temp):
      os.makedirs(self.build_temp)

    subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
    subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)
