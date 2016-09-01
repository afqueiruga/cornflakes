from cornflakes import *
from husk_pair_bond import *

import numpy as np

gdim = 2
Nx = 20
Ny = 20
W = 1.0
H = 1.0
Rad = 0.1
K = 1000.0
L0 = 0.05

# Init the positions
X = PP.init_grid(Nx,Ny,[0.0,0.0],[W,0.0],[0.0,H],0.0)
# Make a graph of bonds
HBond = Hypergraph()
cflib.Build_Pair_Graph(HBond.hg, X, Rad)
# Make more data
v = np.zeros(X.shape,dtype=np.double)
params = np.array([ K, L0 ], dtype=np.double)
# Make a Dofmap object to map vertices to DOFs that indices into x, v, R and K
dmap_ptvec = Dofmap_Strided(gdim)
dmap_params = Dofmap_Strided(2) 

# Now, we want some BCs. Let's look for all of the vertices new the bottom
marked_bot = select_nodes(X, lambda x:x[1]<0.0001)
# And get the indices in R,K they're associated with
bcdofs = dmap_ptvec.Get_List(marked_bot)
bcvals = np.zeros(bcdofs.shape,dtype=np.double)
# And lets grab the top, too
marked_top = select_nodes(X, lambda x:x[1]>H-0.0001)
loaddofs = dmap_ptvec.Get_List(marked_top).reshape(len(marked_top),2)[:,0] # Just the x's

# Lets chop up our graph with a filter operation!
import husk_line_test
dmap_yglobal = Dofmap_Strided(4)
htrue,hfalse = Filter(husk_line_test.kernel_line_intersection, HBond,
                      [dmap_ptvec,dmap_yglobal],
                      [X, np.array([0.0,H/2.0, W/2.0,H/2.0],dtype=np.double)])
print htrue.view(), hfalse.view()
assert( htrue.view()[0].shape[0] + hfalse.view()[0].shape[0] == HBond.view()[0].shape[0] )
HBond = hfalse

# Assemble a matrix and a vector
R,K = Assemble_Targets(kernel_linear_spring,HBond,
                       [dmap_ptvec,dmap_params],[X,v,params],
                       X.size)
# Since its in scipy, you can now just do
# from matplotlib import pylab as plt
# plt.spy(K)
# plt.show()
# And look at your matrix! Do other stuff to it to. Lets go into IPython:
# from IPython import embed
# embed()
# And play around with your matrix!

# Add the point-wise forces
R[loaddofs] -= 200.0
# Lets chop up the matrix to put on BCs
Apply_BC(bcdofs,bcvals,  K,R) # Warning, this is slow!

# Solve the matrix however you want
import scipy.sparse.linalg as splin
u = splin.spsolve(K,R)
u = u.reshape(X.shape) # it comes out flat

# Write it out
GraphIO.write_graph("springs.vtk", HBond, X+u,
                    nodefields=[("x",X),("R",R.reshape(X.shape)),("u",u),("v",v)])

print "The average x-displacement at the top was ", u.flatten()[loaddofs].mean()
