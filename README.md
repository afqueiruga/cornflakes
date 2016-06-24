Cornflakes
==========

2015 Alejandro F Queiruga

Lawrence Berkeley Lab

Intro
-----

Cornflakes is a simple serial library for assembling
vectors and matrices from the kernels created by kerngen.
It's meant for super-quick prototyping for small problems
allowing for interactivity and inspection of the matrices,
etc.

It will have a sister program that also using kerngen but
runs massively-parallel.

Requirements
------------

Kerngen in path, scipy, numpy. Output is to vtk files.

Installing
--------

1. mkdir build
2. cd buid
3. ccmake ..
4. set CMAKE_INSTALL_PREFIX to somewhere local
5. Set USE_LIS or USE_PETSC
6. make
7. make test
8. make install

### OSX Specific:

Things on OSX are tricky. You'll need to know where to link things:

1. set CMAKE_C_COMPILER and CMAKE_CXX_COMPILER gcc-mp-5
2. set SWIG_EXECUTABLE:
/opt/local/bin/swig 
3. Set EXTRA_INCLUDES:
port contents numpy | grep include
4. Set PYTHON_INCLUDE_DIR
port contents python27 | grep include
5. Set PYTHON_LIBRARY
port contents python27 | grep libpython2.7
