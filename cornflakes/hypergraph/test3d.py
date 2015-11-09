#print "1"
from pyHypergraph import Hypergraph
import numpy as np
import scipy.sparse.linalg as splin
import mylibrary
import particle_placers
from dofmap import *
from Assemble import *
import graphio
#print "A"
pks=mylibrary.cvar.particle3d_kernel_strct
print "B"
X = particle_placers.init_cube(11,11,11, [0,0,0], [10,0,0],[0,10,0],[0,0,10], 0.0)
hole = particle_placers.sphere_test(np.array((5.0,5.0,0.0)),10.0)
X = particle_placers.delete_particles(X,[hole])
print "B,5"
hyper = Hypergraph()
mylibrary.Build_Particle_Graph(hyper.hg, X, 1.5)
print "C"
dofmap = make_dofmap(hyper, pks, X.shape[1])
dofnew = DofMap(hyper,pks,X.shape[0])

R = np.zeros(X.size,dtype=np.double)
x = X.copy()
#x[:,0]*=1.1 #*x[:,0]
vel = np.zeros(X.shape)
params = np.array([1.0,-2.0],dtype=np.double)

botnodes = select_nodes(X,lambda a:a[1]<0.1)
topnodes = select_nodes(X,lambda a:a[1]>9.9)
botdofs = np.hstack([dofnew.vertex_dofs[botnodes].flatten(),
                     dofnew.vertex_dofs[topnodes][:,0],
                     dofnew.vertex_dofs[topnodes][:,2]])
botvals = np.zeros(botdofs.shape, dtype=np.double)
topdofs = dofnew.vertex_dofs[topnodes][:,1]
print topdofs
def DYNAMICS():
    from RKnew import RKbase, exRK, imRK
    def sys_mech(time,tang=False):
        #global R
        if tang:
            R = Assemble_Targets((1,),pks,hyper.hg,dofmap,x.size, [x,vel,X,params])[0]
            R[topdofs] += 1.0
            return R,None,None
        else:
            R = Assemble_Targets((1,),pks,hyper.hg,dofmap,x.size, [x,vel,X,params])[0]
            R[topdofs] += 1.0
            return R
    def bcapp_mech(K,R,t,hold=False):
        #embed()
        Apply_BC(botdofs, botvals, K,R)
    def update():
        pass
    Mechfield = RKbase.RK_field(2,[vel.ravel(),x.ravel()],None,
                           sys_mech,bcapp_mech,update)
    Tmax = 10.0
    NT = 1000
    h = Tmax/NT

    step = exRK.exRK(h, exRK.exRK_table['FWEuler'], [Mechfield])
    #step = imRK.DIRK(h, imRK.LDIRK['BWEuler'], [Mechfield])
    from IPython import embed
    #embed()
    print "whu"
    
    print "uh"
    output()
    #step.march()
    NT=1000
    for tx in xrange(NT):
        print "At ",tx
        step.march()
        if tx % 10 ==0:
            output()
outcnt=0
def output():
    global outcnt
    graphio.write_graph("outs/foo_{0}.vtk".format(outcnt),hyper,x, {'x':X,'v':vel})
    outcnt += 1

DYNAMICS()
