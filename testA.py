from cornflakes import *
import numpy as np
import scipy

gdim = 2

#
# Place the particles and build a connectivity
#
hyper = Hypergraph()
X = PP.init_grid(10,10,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.0)
hole = PP.sphere_test(np.array((5.0,5.0)),2.0)
X = PP.delete_particles(X,[hole])
cflib.Build_Pair_Graph(hyper.hg, X, 2.5)

Npart = X.shape[0]
hyper.Add_Edge_Vertex(Npart)
#
# Set up arrays
#
x = X.copy()
v = np.zeros(X.shape)

params = np.array([20.0, -2.0], dtype=np.double)
#params = np.zeros((hyper.hg.he.n_edge,2), dtype=np.double)
#params[:,0] = 20.0
#params[:,1] = -2.0
alpha = np.ones([ hyper.hg.he.n_edge ], dtype=np.double)

# Set up the dofmaps
dmap_vec = Dofmap_Strided(gdim)
dmap_bond = Dofmap_Strided(1,-X.shape[0])
dmap_global = Dofmap_Strided(1)

# The lists
# TODO: FRINGE OF DOF_GLOBAL IN ASSEM
fields = [x,v,X,alpha,params]
dmaps  = [dmap_vec,dmap_bond,dmap_global]

ke = cflib.cvar.kern_peri
#
# Select DOFs for boundary conditions
#
botnodes = select_nodes(X,lambda a:a[1]<0.1)
topnodes = select_nodes(X,lambda a:a[1]>9.9)

from IPython import embed
#embed()

botdofs = dmap_vec.Get_List(botnodes)
botvals = np.zeros(botdofs.shape, dtype=np.double)
topdofs = dmap_vec.Get_List(topnodes)
topdofs = topdofs.reshape((topdofs.size/2,2))[:,1]
topvals = np.zeros(topdofs.shape, dtype=np.double)

embed()


def DYNAMICS():
    import cornflakes.RKnew
    from cornflakes.RKnew import RKbase, exRK, imRK
    def sys_mech(time,tang=False):
        #global R
        if tang:
            R,KX,KV =  Assemble_Targets(ke, hyper, dmaps,fields, X.size)
            #R[topdofs] += 15.0
            #R*=-1.0
            #embed()
            return R,KX,KV
        else:
            R,KX,KV =  Assemble_Targets(ke, hyper, dmaps,fields, X.size)
            #R[topdofs] += 5.0
            #R*=-1.0
            return R
    def bcapp_mech(K,R,t,hold=False):
        Apply_BC(botdofs, botvals, K,R)
        Apply_BC(topdofs, topvals, K,R)
    def update():
        pass
    Mechfield = RKbase.RK_field(2,[v.ravel(),x.ravel()], scipy.sparse.eye(X.size).tocsr(),
                           sys_mech,bcapp_mech,update)
    Tmax = 1.0
    NT = 40
    h = Tmax/NT
    #v.ravel()[topdofs]= 0.01*10.0/Tmax
    v.ravel()[topdofs]= 0.01*10.0/Tmax
    #step = exRK.exRK(h, exRK.exRK_table['RK4'], [Mechfield])
    step = imRK.DIRK(h, imRK.LDIRK['BWEuler'], [Mechfield])
    output()
    for tx in xrange(NT):
        step.march()
        if tx % 1 ==0:
            output()
        #R = Assemble_Targets((1,),damagekernel, hyper.hg,
        #                     dofmap_edges, hyper.hg.n_edge, [x,vel,X,alpha,params])[0]
        #alpha[:] = R[:]
    



outcnt=0
def output():
    global outcnt
    #GraphIO.write_graph("outs/foo_{0}.vtk".format(outcnt),hyper,x, [('x',X),('v',v)],{'alpha':alpha})
    GraphIO.write_silo("out/foo_{0}.silo".format(outcnt),hyper,x, cycle=outcnt,time=outcnt,
                       nodefields=[("x",x),("v",v)], edgefields=[("alpha",alpha)])
    outcnt += 1

DYNAMICS()
#GraphIO.write_silo_meshfile("out/foo_mesh.silo", hyper,X)
#for t in xrange(10):
#    GraphIO.write_silo_datfile("out/foo_{0}.silo".format(t),"foo_mesh.silo",cycle=t,time=t, nodefields=[("m2",m2),("x",v),("v",X),("m",m)], edgefields=[("alpha",alpha)])
#    m2[:]*=2.0
