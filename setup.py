import sys
from distutils.core import setup
from distutils.extension import Extension

USE_CYTHON = False
if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension('_hotcent', ['hotcent/extensions' + ext]),
              ]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, annotate=True)

setup(
  name='Hotcent',
  ext_modules=extensions,
)
