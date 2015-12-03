from cornflakes import *
import numpy as np
import scipy

gdim = 2
Rad = 1.5
ke = cflib.cvar.kern_state

#
# Place the particles and build a connectivity
#

X = PP.init_grid(11,11,[0.0,0.0],[10.0,0.0],[0.0,10.0],0.0)
#hole = PP.sphere_test(np.array((5.0,5.0)),2.0)
#X = PP.delete_particles(X,[hole])

hyper = Hypergraph()
cflib.Build_Adjacency_Graph_Uniform(hyper.hg, X, Rad)

hyperbond = Hypergraph()
cflib.Build_Pair_Graph(hyperbond.hg, X,Rad)

Npart = X.shape[0]
#hyper.Add_Edge_Vertex(Npart)

y = X.copy()
v = np.zeros(X.shape)
params = np.array([Rad, 1.0, -0.5 ,1.0], dtype=np.double)

dmap_vec    = Dofmap_Strided(gdim)
dmap_global = Dofmap_Strided(1)

fields = [ y,v,X, params ]
dmaps  = [ dmap_global, dmap_vec ]

#
# Select DOFs for boundary conditions
#
botnodes = select_nodes(X,lambda a:a[1]<0.1)
topnodes = select_nodes(X,lambda a:a[1]>9.9)

from IPython import embed
#embed()

botdofs = dmap_vec.Get_List(botnodes)
botdofs = botdofs.reshape((botdofs.size/2,2))[:,1]
botvals = np.zeros(botdofs.shape, dtype=np.double)
topdofs = dmap_vec.Get_List(topnodes)
topdofs = topdofs.reshape((topdofs.size/2,2))[:,1]
topvals = np.zeros(topdofs.shape, dtype=np.double)


def DYNAMICS():
    import cornflakes.RKnew
    from cornflakes.RKnew import RKbase, exRK, imRK
    def sys_mech(time,tang=False):
        #global R
        if tang:
            R,KX,KV =  Assemble_Targets(ke, hyper, dmaps,fields, X.size)
            return R,KX,KV
        else:
            R =  Assemble_Targets(ke, hyper, dmaps,fields, X.size)[0]
            #embed()
            return R
    def bcapp_mech(K,R,t,hold=False):
        #embed()
        
        Apply_BC(botdofs, botvals, None,R)
        Apply_BC(topdofs, topvals, None,R)
    def update():
        pass
    Mechfield = RKbase.RK_field(2,[v.ravel(),y.ravel()], None,#scipy.sparse.eye(X.size).tocsr(),
                           sys_mech,bcapp_mech,update)
    Tmax = 15.0
    NT = 1000
    h = Tmax/NT
    #v.ravel()[topdofs]= 0.01*10.0/Tmax
    v.ravel()[topdofs]= 1.0/Tmax
    step = exRK.exRK(h, exRK.exRK_table['RK2-mid'], [Mechfield])
    #step = imRK.DIRK(h, imRK.LDIRK['BWEuler'], [Mechfield])
    output()
    for tx in xrange(NT):
        print tx
        step.march()
        if tx % 10 ==0:
            output()
outcnt=0
def output():
    global outcnt
    #GraphIO.write_graph("outs/foo_{0}.vtk".format(outcnt),hyper,x, [('x',X),('v',v)],{'alpha':alpha})
    GraphIO.write_silo("out/state_{0}.silo".format(outcnt),hyperbond,y, cycle=outcnt,time=outcnt,
                       nodefields=[("x",X),("v",v)], edgefields=[])
    outcnt += 1

DYNAMICS()
