from cornflakes import *
import numpy as np
import scipy
from scipy import sparse as sp
from scipy import linalg
from scipy.sparse import linalg as splin
gdim = 2
#Rad = 1.
ke = cflib.cvar.kernel_darcy_afq
kejac = cflib.cvar.kernel_darcy_afq_CG
kesup = cflib.cvar.kernel_darcy_support_afq
dom_w = 2.0
dom_h = 2.0

#
# Place the particles and build a connectivity
#
outcnt = 0
def Do_Sim(Rad,Nside, pold=None, Xold=None):
    # Build particles and connectivity
    X = PP.init_grid(Nside+int(Nside*(2.0*Rad/dom_h)),Nside+int(Nside*(2.0*Rad/dom_h)),
                     [-dom_w/2.0 - Rad,-dom_h/2.0 - Rad],
                     [ dom_w + 2.0*Rad,0.0],
                     [0.0,dom_h+2.0*Rad],0.0)
    #hole = PP.sphere_test(np.array((5.0,5.0)),2.0)
    #X = PP.delete_particles(X,[hole])
    hyper = Hypergraph()
    cflib.Build_Pair_Graph(hyper.hg, X,1.1*Rad)
    Npart = X.shape[0]

    area =  dom_w/(Nside-1.0)*dom_h/(Nside-1.0)
    
    # Initial Dmaps
    dmap_vec    = Dofmap_Strided(gdim)
    dmap_sca    = Dofmap_Strided(1)
    dmap_global = Dofmap_Strided(gdim)
    dmap_param  = Dofmap_Strided(4)
    
    dmaps  = [ dmap_sca, dmap_vec, dmap_param ]

    # Initialize fields
    p = np.zeros(Npart,dtype=np.double)
    Dp = np.ones(Npart,dtype=np.double)
    #params = np.array([Rad, 1.0 ], dtype=np.double)
    params = np.zeros((Npart,4))
    params[:,0]=Rad
    params[:,1]=1.0
    params[:,2] = area

    fields = [ p, Dp, X, params ]


    SUP, = Assemble_Targets(kesup, hyper, dmaps,[p,X,params], Npart)
    params[:,3] = area*(1.0-(area/4.0*1.0/3.0/Rad))**2.0 + SUP[:]
    
    # Make the initial guess if one was provided
    from IPython import embed
    #embed()\
    print " Generating guess"
    if pold!=None and Xold!=None:
        pold2 = pold.view()
        pold2.shape = (pold.shape[0],1)
        p2 = p.view()
        p2.shape = (p.shape[0],1)
        #embed()
        cflib.Interpolate_np(pold2,Xold,
                             p2,X, 1.1*abs( Xold[0,0]-Xold[1,0]) )
        #embed()
    print " Setting up"

    
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

    allnodes = np.hstack( [botnodes, topnodes, rightnodes, leftnodes] )
    alldofs = dmap_sca.Get_List(allnodes)
    allvals = np.zeros(alldofs.shape, dtype=np.double)

    #blockdofs = dmap_sca.Get_List(blocknodes) This was a kludge!
    #params[blocknodes,1]=0.005

    #outcnt = 0
    def output():
        global outcnt
        GraphIO.write_silo("out/darcy_{0}.silo".format(outcnt),hyper,X, cycle=outcnt,time=outcnt,
                       nodefields=[("x",X),("p",p),("Sup",params[:,3])], edgefields=[])
        outcnt += 1

    eps = 1.0
    TOL = 1.0e-10
    itr = 0
    p[alldofs] = 0.0
    def sys(DP):
        R,DR = Assemble_Targets(kejac, hyper, dmaps,[p,DP, X,params],Npart)
        #R[:] += 1.0#sourcenodes]+=1.0
        #R[drainnodes]-=1.0
        Apply_BC(alldofs, allvals, None, DR)
        return DR
    import inspect
    def callback(P):
        frame = inspect.currentframe().f_back
        print(frame.f_locals['resid'])
        pass #p=P
    R,DR =  Assemble_Targets(kejac, hyper, dmaps,
                             [np.zeros(p.shape),np.zeros(p.shape),X,params],
                             Npart)
    R[:] += 1.0#sourcenodes]+=1.0
    R[drainnodes]-=1.0
    Apply_BC(alldofs, allvals, None, R)


    L = splin.LinearOperator( (Npart,Npart), matvec=sys, dtype=np.double)
    print " Solving"
    x,info = splin.cg(L,R,x0=p,callback = None)
    #from IPython import embed
    #embed()
    p = x
    #while itr < 10 and eps > TOL:
    #    R,K =  Assemble_Targets(ke, hyper, dmaps,fields, Npart)
        #from IPython import embed
        #embed()
    #    R[:] += 1.0#sourcenodes]+=1.0
    #    R[drainnodes]-=1.0
    #    Apply_BC(topdofs, topvals, K,R)
    #    Apply_BC(botdofs, botvals, K,R)
    #    Apply_BC(leftdofs, leftvals, K,R)
    #    Apply_BC(rightdofs, rightvals, K,R)
        
    #    deltap = splin.spsolve(K,R)
    #    p -= deltap
    #    eps = linalg.norm(deltap)
    #    print itr, ":",eps
    #    itr += 1
    output()
    return p,X

def error(p,X):
    tot = 0.0
    k=1.0
    print " Calculating error"
    panal = np.zeros(p.shape)
    t = lambda k,x,y : np.sin(k*np.pi*(1.0+x)/2.0)/( k*k * np.sinh( k*np.pi )) * \
        ( np.sinh( k*np.pi*(1.0+y)/2.0 ) + np.sinh( k*np.pi*(1.0-y)/2.0) )
    for i,y in enumerate(X):
        
        panal[i] = (1.0-y[0]**2.0)/2.0 - 16.0/(np.pi)**3 * (
            t(1,y[0],y[1]) + t(3,y[0],y[1]) + t(7,y[0],y[1]) + t(9,y[0],y[1])
        )
    
    return np.linalg.norm(p-panal)*dom_w/np.sqrt(p.shape[0])

pold = None
Xold = None
f_results = open("convergence","a")
for i in [2**j+1 for j in range(3,9) ]:
    for RF in [1.5,2.5,3.5,4.5]:
        print RF," " , i
        r = RF*dom_w/(float(i)-1)
        p,X = Do_Sim(r, i, pold,Xold)
        pold = p
        Xold = X
        e = error(p,X)
        line =  " ".join([str(a) for a in (RF, dom_w/(float(i)-1), np.max(p), np.mean(p),e )])
        print line
        f_results.write( line + "\n" )
        f_results.flush()
#        print np.mean(p)

f_results.close()
