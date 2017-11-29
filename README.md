Cornflakes
==========

2015 Alejandro F Queiruga

Lawrence Berkeley Lab

Intro
-----

Cornflakes is a simple library for assembling
vectors and matrices from the kernels created by the accompanying
Domain Specific Language popcorn. In a sense, cornflakes
is the runtime to the popcorn DSL. 
It's meant for super-quick prototyping for small problems
allowing for interactivity and inspection of the matrices,
etc.

It is designed to be run in parallel, but, like many things,
I haven't gotten to that part yet. I might just rewrite it 
all in Julia first.

Philosophy
----------

1. The user must retain control of main.
2. Agnostic to how the kernel was written. It just has to have 
the kernel struct. The first DSL was already replaced (kerngen).
3. The code generation should run pre-compile time, not just-in-time.
It gets expensive. Some of my kernels take 5 mintues to run. The 
production code should be given already compiled shared object files.

Example
-------


Requirements
------------

Popcorn in path, scipy, numpy. Output is to vtk files.

