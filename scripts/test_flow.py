from cornflakes import *
import numpy as np
import scipy
from scipy import sparse as sp
from scipy import linalg
from scipy.sparse import linalg as splin
gdim = 2
#Rad = 1.
ke = cflib.cvar.kernel_darcy_afq
dom_w = 2.0
dom_h = 2.0

#
# Place the particles and build a connectivity
#
outcnt = 0
def Do_Sim(Rad,Nside):
    # Build particles and connectivity
    X = PP.init_grid(Nside,Nside,
                     [-dom_w/2.0,-dom_h/2.0],
                     [ dom_w,0.0],
                     [0.0,dom_h],0.0)
    #hole = PP.sphere_test(np.array((5.0,5.0)),2.0)
    #X = PP.delete_particles(X,[hole])
    hyper = Hypergraph()
    cflib.Build_Pair_Graph(hyper.hg, X,Rad)
    Npart = X.shape[0]

    # Initial Dmaps
    dmap_vec    = Dofmap_Strided(gdim)
    dmap_sca    = Dofmap_Strided(1)
    dmap_global = Dofmap_Strided(gdim)
    dmaps  = [ dmap_sca, dmap_vec, dmap_vec ]
    
    # Initialize fields
    p = np.zeros(Npart,dtype=np.double)
    #params = np.array([Rad, 1.0 ], dtype=np.double)
    params = np.zeros(X.shape)
    params[:,0]=Rad
    params[:,1]=1.0
    fields = [ p, X, params ]

    # Select DOFs for boundary conditions
    botnodes   = select_nodes(X,lambda a:a[1]<-dom_w/2.0 + dom_w/(Nside+1))
    topnodes   = select_nodes(X,lambda a:a[1]> dom_w/2.0 - dom_w/(Nside+1))
    rightnodes = select_nodes(X,lambda a:a[0]<-dom_h/2.0 + dom_h/(Nside+1))
    leftnodes  = select_nodes(X,lambda a:a[0]> dom_h/2.0 - dom_h/(Nside+1))
#    topnodes = select_nodes(X,lambda a:a[1]>10.0- 10.0/(Nside+1)) # and a[0]>5.0-Rad and a[0]<5.0+Rad)
    blocknodes = select_nodes(X, lambda a: (a[1]>1 and a[1]<9) and (a[0]>3 and a[0]<7) )
    #sourcenodes = select_nodes(X, lambda a: (a[0]**2+a[1]**2)<1.0 )
    #drainnodes = select_nodes(X, lambda a: ((10.0-a[0])**2+(10.0-a[1])**2)<1.0 )
    sourcenodes = select_nodes(X, PP.sphere_test(np.array((2.0,5.0)),1.0) )
    drainnodes = []# select_nodes(X, PP.sphere_test(np.array((8.0,5.0)),1.0) )

    botdofs = dmap_sca.Get_List(botnodes)
    botvals = np.zeros(botdofs.shape, dtype=np.double)
    topdofs = dmap_sca.Get_List(topnodes)
    topvals = np.zeros(topdofs.shape, dtype=np.double)
    rightdofs = dmap_sca.Get_List(rightnodes)
    rightvals = np.zeros(rightdofs.shape, dtype=np.double)
    leftdofs = dmap_sca.Get_List(leftnodes)
    leftvals = np.zeros(leftdofs.shape, dtype=np.double)
    #blockdofs = dmap_sca.Get_List(blocknodes) This was a kludge!
    #params[blocknodes,1]=0.005

    #outcnt = 0
    def output():
        global outcnt
        GraphIO.write_silo("out/darcy_{0}.silo".format(outcnt),hyper,X, cycle=outcnt,time=outcnt,
                       nodefields=[("x",X),("p",p)], edgefields=[])
        outcnt += 1

    eps = 1.0
    TOL = 1.0e-10
    itr = 0
    p[topdofs] = 0.0
    while itr < 10 and eps > TOL:
        R,K =  Assemble_Targets(ke, hyper, dmaps,fields, Npart)
        #from IPython import embed
        #embed()
        R[:] += 1.0#sourcenodes]+=1.0
        R[drainnodes]-=1.0
        Apply_BC(topdofs, topvals, K,R)
        Apply_BC(botdofs, botvals, K,R)
        Apply_BC(leftdofs, leftvals, K,R)
        Apply_BC(rightdofs, rightvals, K,R)
        
        deltap = splin.spsolve(K,R)
        p -= deltap
        eps = linalg.norm(deltap)
        print itr, ":",eps
        itr += 1
    output()
    return p,X

def error(p,X):
    tot = 0.0
    k=1.0
    panal = (10.0*X[:,1]-X[:,1]*X[:,1])/2.0 #((50+k)*X[:,1]-5.0*X[:,1]*X[:,1])/(10.0*k)
    
    return np.linalg.norm(p-panal)
for i in [13,24,43,83,107]:
    r = 4.0*dom_w/(float(i)-1)
    p,X = Do_Sim(r, i)
#    print error(p,X)
    print np.mean(p)
    

