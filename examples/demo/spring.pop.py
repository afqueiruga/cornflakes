from popcorn import *
# This is an example of a spring-law, the simplest possible kernel
#
# P0 ---- P1
#
gdim = 2
# DofMap Type-defs
PtSca = DofSpace(   1,0,2) # a sclar field associated with P0 and P1 ( 2 numbers )
PtVec = DofSpace(gdim,0,2) # a vector field associated with P0 and P1 ( 4 numbers )
Param = DofSpace(2,-1) # a field to input global parameters. 2 numbers, not associated with a vertex
# Inputs and outputs
i_x = Input("x", PtVec)
i_params = Input("params",Param)
o_R = Output("R", [PtVec], 1)
o_K = Output("K", [PtVec], 2)
# Grabbing symbols
x0,x1 = i_x.Vertex_Handles(0,1) # Split into blocks of gdim
K, L0 = i_params.Entry_Handles(0,1) # split into each entry in the vector
# Vector algebra
r = x1 - x0
rabs = sqrt((r.T*r)[0,0])
f_expr = K*(rabs-L0) * r/rabs
# Expressions for the contributions
R_expr = Matrix([ f_expr, -f_expr])
K_expr = R_expr.jacobian( i_x.as_matrix() )
# The kernel prgm and object
prgm = [
    Asgn(o_R, R_expr, "+="),
    Asgn(o_K, K_expr, "+=")
    ]
kernel_linear_spring = Kernel("linear_spring",
                              [i_x, i_params],
                              [o_R, o_K],
                              listing=prgm)
# Write a husk
Husk("pair_bond", [ kernel_linear_spring ])
