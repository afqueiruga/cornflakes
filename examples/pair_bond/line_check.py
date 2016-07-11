from kerngen import *
# A routine to see if two lines are touching
# [P1 P2 C1]
# y is a global

gdim = 2
xPtVec = DofSpace(gdim, 0,2)
xCellSca = DofSpace(1, 2,3)
yVec = DofSpace(2*gdim, -1)


x_d = Input("x",xPtVec)
y_d = Input("y",yVec)
y_d.dim = 2*gdim # TODO: Fix in kerngen

x1,x2 = x_d.Vertex_Split()
y1,y2 = Matrix(y_d.Handle()[0:2]), Matrix(y_d.Handle()[2:4])

s,r = symbols('s r')
itscn = solve( x1+ (x2-x1)*s - (y1+(y2-y1)*r) , s,r)

Testout = TensorVariable("testout",1,1)
prgm = [
    IfStmt( And(Ge(itscn[s],0),
                Ge(itscn[r],0),
                Le(itscn[s]**2, ((x2-x1).T*(x2-x1))[0,0] ),
                Le(itscn[r]**2, ((y2-y1).T*(y2-y1))[0,0] ) ),
            Assignment(Testout, Matrix([[1]]) )
            )
    ]
outp_Testout = Output("testout",[xCellSca],1)

kernel_testout = Kernel("line_intersection",
                        [x_d,y_d],[outp_Testout],
                        listing=prgm)

modname = "line_test"
pylink(modname,
       [ kernel_testout ],
       targetdir="husk_"+modname,config=boilerplates.pylink.config_osx)
