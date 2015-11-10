import numpy as np
import scipy.sparse.linalg as splin

from pyHypergraph import Hypergraph
import mylibrary
import particle_placers
from dofmap import *
import graphio
from Assemble import *

#
# Place the particles and build a connectivity
#
hyper = Hypergraph()
X = particle_placers.init_grid(50,50,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.1)
hole = particle_placers.sphere_test(np.array((5.0,5.0)),2.0)
X = particle_placers.delete_particles(X,[hole])
mylibrary.Build_Particle_Graph(hyper.hg, X, 0.4)

#
# Handles to the kernels
#
pks=mylibrary.cvar.particle_kernel_strct
damagekernel=mylibrary.cvar.damage_kernel_strct


#
# Make DOF mappings
#
dofmap = make_dofmap(hyper, pks, X.shape[1])
dofnew = DofMap(hyper,pks,X.shape[0])
dofmap_edges = np.arange(hyper.hg.n_edge,dtype=np.intc)

#
# Select DOFs for boundary conditions
#
botnodes = select_nodes(X,lambda a:a[1]<0.1)
topnodes = select_nodes(X,lambda a:a[1]>9.9)
#botdofs = np.hstack([dofnew.vertex_dofs[botnodes].flatten(),
#                     dofnew.vertex_dofs[topnodes][:,0]])
botdofs = np.hstack([dofnew.vertex_dofs[botnodes][:,1],
                     dofnew.vertex_dofs[topnodes][:,1]])
botvals = np.zeros(botdofs.shape, dtype=np.double)
topdofs = dofnew.vertex_dofs[topnodes][:,1]
topvals = np.zeros(topdofs.shape, dtype=np.double)


#
# Initialize the fields
#
x = X.copy()
vel = np.zeros(X.shape)
params = np.array([20.0,-2.0],dtype=np.double)
alpha = np.ones((hyper.hg.n_edge,),dtype=np.double)
#x[:,0]*=1.1 #*x[:,0]


def DYNAMICS():
    from RKnew import RKbase, exRK, imRK
    def sys_mech(time,tang=False):
        #global R
        if tang:
            R,KX,KV = Assemble_Targets((1,2,2), pks,hyper.hg,dofmap,x.size, [x,vel,X,alpha,params])
            #R[topdofs] += 5.0
            #R*=-1.0
            return R,KX,KV
        else:
            R,KX,KV = Assemble_Targets((1,2,2), pks,hyper.hg,dofmap,x.size, [x,vel,X,alpha,params])
            #R[topdofs] += 5.0
            #R*=-1.0
            return R
    def bcapp_mech(K,R,t,hold=False):
        Apply_BC(botdofs, botvals, K,R)
        Apply_BC(topdofs, topvals, K,R)
    def update():
        pass
    Mechfield = RKbase.RK_field(2,[vel.ravel(),x.ravel()], scipy.sparse.eye(X.size).tocsr(),
                           sys_mech,bcapp_mech,update)
    Tmax = 100.0
    NT = 40
    h = Tmax/NT
    vel.ravel()[topdofs]= 0.01*10.0/Tmax
    
    #step = exRK.exRK(h, exRK.exRK_table['RK4'], [Mechfield])
    step = imRK.DIRK(h, imRK.LDIRK['BWEuler'], [Mechfield])
    output()
    for tx in xrange(NT):
        step.march()
        if tx % 1 ==0:
            output()
        R = Assemble_Targets((1,),damagekernel, hyper.hg,
                             dofmap_edges, hyper.hg.n_edge, [x,vel,X,alpha,params])[0]
        alpha[:] = R[:]
    
outcnt=0
def output():
    global outcnt
    graphio.write_graph("outs/foo_{0}.vtk".format(outcnt),hyper,x, {'x':X,'v':vel},{'alpha':alpha})
    outcnt += 1

DYNAMICS()
