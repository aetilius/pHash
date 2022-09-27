import os
from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

phash_ext = Extension('phash',
                      sources=['phash.pyx', os.path.abspath('../../src/pHash.cpp')],
                      libraries=['png', 'tiff'],
                      language='c++')

setup(
    ext_modules=cythonize([phash_ext], language_level='3'))
