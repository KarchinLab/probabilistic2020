from distutils.core import setup
from distutils.extension import Extension
import sys

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, 'permutation2020/cython'))

SRC_DIR = 'permutation2020'

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.cpp'
extensions = [
    Extension(SRC_DIR + ".cython.uniform_kde",
              [SRC_DIR +'/cython/uniform_kde'+ext],
              language='c++',
              include_dirs=[SRC_DIR + '/cython/']),
                            # np.get_include()]),
    Extension(SRC_DIR + ".cython.gaussian_kde",
              [SRC_DIR + '/cython/gaussian_kde'+ext],
              language='c++',
              include_dirs=[SRC_DIR + '/cython/']),
                            # np.get_include()]),
    Extension(SRC_DIR + ".cython.cutils",
              [SRC_DIR + "/cython/cutils"+ext],
              language='c++',
              include_dirs=[SRC_DIR + '/cpp/',
                            SRC_DIR + '/cython/'])
                            #np.get_include()])
]


if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

if 'build_ext' in sys.argv:
    # just build cython extension module if build_ext subcommand is used
    setup(ext_modules = extensions)
else:
    import permutation2020
    version = permutation2020.__version__
    AUTHOR = 'Collin Tokheim'
    EMAIL = 'fake@gmail.com'
    URL = 'https://github.com/ctokheim/2020permutation'
    DESCRIPTION = '20/20 Permutation Test'
    PACKAGES = [SRC_DIR, SRC_DIR + '.python', SRC_DIR + '.cython', SRC_DIR + '.cpp']
    setup(name='permutation2020',
          version=version,
          description=DESCRIPTION,
          author=AUTHOR,
          author_email=EMAIL,
          url=URL,
          packages=PACKAGES,
          install_requires=['numpy', 'scipy', 'pandas==0.12.0', 'pysam'],
          scripts=['bin/permutation_test.py', 'bin/extract_gene_seq.py'],
          long_description=open('README.rst').read(),
          classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                       'Environment :: Console',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research'],
          ext_modules=extensions
    )
