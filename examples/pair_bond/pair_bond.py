from kerngen import *
# This is an example of a spring-law, the simplest possible kernel
# The edge is fixed-length and looks like this:
#
# P0 ---- P1
#

#
# The inputs set up
#
# This will generate any dimension
gdim = 2

# Set up the DofSpaces. These say what sort of data lives on the edge.
# These are sort-of Typedefs
PtSca = DofSpace(   1,0,2) # a sclar field associated with P0 and P1 ( 2 numbers )
PtVec = DofSpace(gdim,0,2) # a vector field associated with P0 and P1 ( 4 numbers )
Global = DofSpace(2,-1) # a field to input global parameters. 2 numbers, not associated with a vertex

# Now register inputs
x_d = Input("x", PtVec)
v_d = Input("v", PtVec)
params_d = Input("params",Global)

# So those are kerngen types, we need SymPy types. There are a few ways to grab them.
# Get a matrix of all of the positions
x_N = x_d.Handle()#.reshape(2,2)
# Get a list of matrices of each entry (and then tuple-unpack it)
x0, x1 = x_d.Vertex_Split()
v0, v1 = v_d.Vertex_Split()
# Get just one entry in the array. (Not by Vertex_Number, but by array index! Always a scalar!)
K = params_d.Entry_Handle(0) # Spring constant
L0 = params_d.Entry_Handle(1) # Resting length (all bonds the same length in this example)

#
# Now define some maths. This is all in the Sympy Matrix class
#
r = x1 - x0
rabs = sqrt((r.T*r)[0,0]) # Ew! Sympy doesn't treat 1x1 matrices nicely,
# so sometimes you have to do this to grab a scalar-sympy expression
# E.g., when you do a dot-product

# The force law
f_expr = K*(rabs-L0) * r/rabs
# Now just stack 'em. The ouput looks like
# 
#      fx
# R =  fy
#     -fx
#     -fy
#
R_expr = Matrix([ f_expr, -f_expr ])
K_expr = R_expr.jacobian(x_N)

# Now we have to write the "program" that does the calculation.
# It seems silly for this example, but this feature is what lets you write
# kernels that have complicating stuff, like an iterative method, or variable
# length inputs
R = TensorVariable("R",1,4) # The "Load" vector. rank-1 (vector), 4 entries
K = TensorVariable("K",2,4) # The "Tangent Matrix. rank-2, (matrix), 4x4 entries
prgm = [
    Assignment(R, R_expr),
    Assignment(K, K_expr)
    ]

#
# The final kernel definition
#
R_outp = Output("R", [PtVec], 1) # Rank one, a vector associated with the points.
K_outp = Output("K", [PtVec], 2)
# The string name has to match with the TensorVariable. Yes, VERY UGLY! I WILL FIX IT
# The difference: TensorVariables are intermediate variables, e.g. you can declare
# a vector or\ matrix of tempory data that you operate on iteratively or something
# that isn't an input or an output. Then, you need to Output type to say how the
# kernel's output needs to be placed into global Matrix/Vectors during assembly
# by cornflakes.

kernel_linear_spring = Kernel("linear_spring", # The name. 
                              [x_d,v_d,params_d], # List of inputs
                              [R_outp, K_outp], # List of outputs
                              listing=prgm) # The prgm above

# Write a husk
modname = "husk_pair_bond" # This will name the folder
pylink("pair_bond",
       [ kernel_linear_spring ],
       targetdir=modname, config=boilerplates.pylink.config_osx)
# There's now a folder filled with files. Look at 'em.
# You have the following options:
# 1) cd into the folder and type "make". That'll use swig
#    to turn that folder into a python module.
# 2) Copy and paste the .c and .h files whereever you want
# 3) In your CMakeLists.txt, add_subdirectory(husk_pair_bond) to populate
#    KERNEL_FILES and KERNEL_INCLUDES in your build system
