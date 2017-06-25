def get_include(*args,**kwargs):
  import os
  try:
    from pip import locations
    return os.path.dirname(locations.distutils_scheme('cppmat',*args,**kwargs)['headers'])
  except ImportError:
    return 'include'
