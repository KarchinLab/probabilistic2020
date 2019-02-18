#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
import sys

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, 'prob2020/cython'))

SRC_DIR = 'prob2020'

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False

import numpy as np
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [
    Extension(SRC_DIR + ".cython.cutils",
              [SRC_DIR + "/cython/cutils"+ext],
              language='c++',
              include_dirs=[SRC_DIR + '/cpp/',
                            SRC_DIR + '/cython/',
                            np.get_include()])
]


if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

if 'build_ext' in sys.argv:
    # just build cython extension module if build_ext subcommand is used
    setup(ext_modules = extensions)
else:
    import prob2020
    version = prob2020.__version__
    AUTHOR = 'Collin Tokheim'
    EMAIL = 'fake@gmail.com'
    URL = 'https://github.com/KarchinLab/probabilistic2020'
    DESCRIPTION = 'Probabilistic 20/20'
    PACKAGES = [SRC_DIR, SRC_DIR + '.python',
                SRC_DIR + '.cython', SRC_DIR + '.cpp',
                SRC_DIR + '.console']
    setup(name='probabilistic2020',
          version=version,
          description=DESCRIPTION,
          author=AUTHOR,
          author_email=EMAIL,
          url=URL,
          packages=PACKAGES,
          license='Apache License, Version 2.0',
          install_requires=['numpy', 'scipy', 'pandas', 'pysam'],
          package_data={
              SRC_DIR+'.console': ['*.R']
          },
          entry_points={
              'console_scripts':[
                  'probabilistic2020 = prob2020.console.probabilistic2020:cli_main',
                  'mut_annotate = prob2020.console.annotate:cli_main',
                  'extract_gene_seq = prob2020.console.extract_gene_seq:cli_main',
                  'simulate_non_silent_ratio = prob2020.console.simulate_non_silent_ratio:cli_main'
              ]
          },
          long_description=open('README.rst').read(),
          classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                       'Environment :: Console',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research'],
          ext_modules=extensions
    )
