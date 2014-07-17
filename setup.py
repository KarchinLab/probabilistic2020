from distutils.core import setup
from distutils.extension import Extension
import sys
# from Cython.Distutils import build_ext

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, 'src/cython'))

import numpy as np

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
extensions = [
    #Extension("src.cython.typedefs",
              #['src/cython/typedefs'+ext],
              #language='c++',
              #include_dirs=['src/cython/',
                            #np.get_include()]),
    Extension("src.cython.uniform_kde",
              ['src/cython/uniform_kde'+ext],
              language='c++',
              include_dirs=['src/cython/',
                            np.get_include()]),
    Extension("src.cython.gaussian_kde",
              ['src/cython/gaussian_kde'+ext],
              language='c++',
              include_dirs=['src/cython/',
                            np.get_include()]),
    Extension("src.cython.cutils",
              ["src/cython/cutils"+ext],
              language='c++',
              include_dirs=['src/cpp/', 'src/cython/',
                            np.get_include()])
]


if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

for extension in extensions:
    setup(
        ext_modules = [extension]
    )
