from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = ROOT

ext = Extension(
    name="itrigamma",
    sources=[
        os.path.join(ROOT, "citrigamma.pyx"),
        os.path.join(ROOT, "citrigamma.c"),
    ],
    include_dirs=[np.get_include(), ROOT],
    extra_compile_args=["-O3"],
)

setup(
    name="itrigamma",
    version="0.0.1a",
    ext_modules=cythonize([ext], compiler_directives={"language_level": "3"}),
)