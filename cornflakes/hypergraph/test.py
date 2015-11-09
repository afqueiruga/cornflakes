print "1"

import numpy as np
import scipy.sparse.linalg as splin

print "A"
from pyHypergraph import Hypergraph
h = Hypergraph(2)
h.Push_Edge(np.array([1,2],dtype=np.intc))
print "B"
import mylibrary
import particle_placers
hyper = Hypergraph()
X = particle_placers.init_grid(21,21,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.0)
hole = particle_placers.sphere_test(np.array((5.0,5.0)),2.0)
X = particle_placers.delete_particles(X,[hole])
print "C"
mylibrary.Build_Particle_Graph(hyper.hg, X, 1.5)
print "D"

pks=mylibrary.cvar.particle_kernel_strct

from dofmap import * 

dofmap = make_dofmap(hyper, pks, X.shape[1])
dofnew = DofMap(hyper,pks,X.shape[0])

R = np.zeros(X.size,dtype=np.double)
x = X.copy()
#x[:,0]*=1.1 #*x[:,0]
vel = np.zeros(X.shape)



params = np.array([3.0],dtype=np.double)
import graphio
from Assemble import *


botnodes = select_nodes(X,lambda a:a[1]<0.1)
topnodes = select_nodes(X,lambda a:a[1]>9.9)
botdofs = np.hstack([dofnew.vertex_dofs[botnodes].flatten(),dofnew.vertex_dofs[topnodes][:,0]])
botvals = np.zeros(botdofs.shape, dtype=np.double)
topdofs = dofnew.vertex_dofs[topnodes][:,1]

def STATICTEST():
    tol = 1e-10
    eps = tol+1.0
    maxiter=10
    itcnt = 0
    output()
    while eps>tol and itcnt < maxiter:
        print "Solving..."
        R,K = Assemble_Vector_Matrix(pks,hyper.hg, dofmap,x.size, [x,vel,X,params])
        Apply_BC(botdofs,botvals, K,R)
        dx = splin.spsolve(K,R)
        eps = np.linalg.norm(dx)
        x.ravel()[:] -= dx[:]
        print "Err = ", eps
        output()

def DYNAMICTEST():
    from RKnew import RKbase, exRK, imRK
    def sys_mech(time,tang=False):
        #global R
        if tang:
            R,KX,KV = Assemble_Targets(pks,hyper.hg,dofmap,x.size, [x,vel,X,params])
            R[topdofs] += 5.0
            #R*=-1.0
            return R,KX,KV
        else:
            R,KX,KV = Assemble_Targets(pks,hyper.hg,dofmap,x.size, [x,vel,X,params])
            #R*=-1.0
            return R
    def bcapp_mech(K,R,t,hold=False):
        #embed()
        Apply_BC(botdofs, botvals, K,R)
    def update():
        pass
    Mechfield = RKbase.RK_field(2,[vel.ravel(),x.ravel()],None,
                           sys_mech,bcapp_mech,update)
    Tmax = 100.0
    NT = 100
    h = Tmax/NT

    #step = exRK.exRK(h, exRK.exRK_table['RK4'], [Mechfield])
    step = imRK.DIRK(h, imRK.LDIRK['BWEuler'], [Mechfield])
    output()
    for tx in xrange(NT):
        step.march()
        if tx % 1 ==0:
            output()
        
outcnt=0
def output():
    global outcnt
    graphio.write_graph("outs/foo_{0}.vtk".format(outcnt),hyper,x, {'x':X,'v':vel})
    outcnt += 1

#R,KX,KV=Assemble_Targets(pks,hyper.hg, dofmap,x.size, [x,vel,X,params])


DYNAMICTEST()
#STATICTEST()
