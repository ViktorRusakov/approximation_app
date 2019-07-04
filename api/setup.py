from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("algebra_grading.pyx", annotate=True),
)

# python setup.py build_ext --inplace