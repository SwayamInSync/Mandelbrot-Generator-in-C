from setuptools import setup, Extension

module = Extension(
    'mandelbrot',
    sources=['main.cpp'],
    libraries=['sleef', 'sleefquad', 'm'],
    extra_compile_args=[
        '-fopenmp', 
        '-Wall',
        '-O3'
    ],
    extra_link_args=['-fopenmp'],
    define_macros=[('_OPENMP', None)]
)

setup(name='mandelbrot',
      version='1.0',
      description='Mandelbrot set calculator',
      ext_modules=[module])
