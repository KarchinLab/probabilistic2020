from distutils.core import setup
from distutils.extension import Extension
import sys
# from Cython.Distutils import build_ext

#setup(
    #cmdclass = {'build_ext': build_ext},
    #ext_modules = [Extension("src.python.cython_utils",
                             #["src/python/cython_utils.pyx"],
                             #language='c++',
                             #include_dirs=['src/cpp/'])]
#)

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [Extension("src.python.cython_utils",
                        ["src/python/cython_utils"+ext],
                        language='c++',
                        include_dirs=['src/cpp/'])]


if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(
    ext_modules = extensions
)
