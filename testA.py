from cornflakes import *
import numpy as np


gdim = 2

#
# Place the particles and build a connectivity
#
hyper = Hypergraph()
X = PP.init_grid(50,50,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.05)
hole = PP.sphere_test(np.array((5.0,5.0)),2.0)
X = PP.delete_particles(X,[hole])
cflib.Build_Pair_Graph(hyper.hg, X, 0.45)

Npart = X.shape[0]
hyper.Add_Edge_Vertex(Npart)
#
# Set up arrays
#
x = X.copy()
v = np.zeros(X.shape)
m = np.arange(X.shape[0],dtype=np.double)
m2 = 2.1*m
#params = np.array([20.0, -2.0], dtype=np.double)
params = np.zeros((hyper.hg.he.n_edge,2), dtype=np.double)
params[:,0] = 20.0
params[:,1] = -2.0
alpha = np.ones([ hyper.hg.he.n_edge ], dtype=np.double)

# Set up the dofmaps
dmap_vec = Dofmap_Strided(gdim)
dmap_bond = Dofmap_Strided(1,-X.shape[0])
dmap_global = Dofmap_Strided(1)

# The lists
# TODO: FRINGE OF DOF_GLOBAL IN ASSEM
fields = [x,v,X,alpha,params]
dmaps  = [dmap_vec,dmap_vec,dmap_vec,dmap_bond,dmap_bond]

ke = cflib.cvar.kern_peri


GraphIO.write_silo_meshfile("out/foo_mesh.silo", hyper,X)
for t in xrange(10):
    GraphIO.write_silo_datfile("out/foo_{0}.silo".format(t),"foo_mesh.silo",cycle=t,time=t, nodefields=[("m2",m2),("x",v),("v",X),("m",m)], edgefields=[("alpha",alpha)])
    m2[:]*=2.0

#R,KX,KV = Assemble_Targets(ke, hyper, dmaps,fields, X.size)


