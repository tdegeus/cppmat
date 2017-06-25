def has_flag(compiler, flagname):
  import tempfile
  with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
    f.write('int main (int argc, char **argv) { return 0; }')
    try:
      compiler.compile([f.name], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
      return False
  return True

def cpp_flag(compiler):
  if   has_flag(compiler,'-std=c++14'): return '-std=c++14'
  elif has_flag(compiler,'-std=c++11'): return '-std=c++11'
  raise RuntimeError('Unsupported compiler: at least C++11 support is needed')

def get_include(*args,**kwargs):
  import os
  try:
    from pip import locations
    return os.path.dirname(locations.distutils_scheme('cppmat',*args,**kwargs)['headers'])
  except ImportError:
    return 'include'
