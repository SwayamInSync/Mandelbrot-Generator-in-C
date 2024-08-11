from setuptools import setup, Extension
import os

module = Extension(
    'mandelbrot',
    sources=['main.cpp'],
    libraries=['sleef', 'sleefquad', 'm', 'tlfloat'],
    include_dirs=[
        '/usr/local/include',
        'install/include'
    ],
    library_dirs=[
        '/usr/local/lib',  # Add the path where 'libtlfloat.dylib' is located
        'install/lib'
    ],
    extra_compile_args=[
        '-fopenmp',
        '-Wall',
        '-O3',
    ],
    extra_link_args=[
        '-fopenmp',
        '-ltlfloat',
        '-lm',
        '-lsleef',
        '-lsleefquad',
        '-Wl,-rpath,/usr/local/lib'
    ]

)

setup(
    name='mandelbrot',
    version='1.0',
    description='Mandelbrot set calculator',
    ext_modules=[module]
)
