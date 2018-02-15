import numpy as np

#
# Generation a Quadrature background mesh
#
Quad_Tables = {
    1: [ (2.0, 0.0) ],
    2: [ (1.0,-np.sqrt(1.0/3.0)),(1.0,np.sqrt(1.0/3.0) ) ],
    3: [ (5.0/9.0, -np.sqrt(3.0/5.0)), (8.0/9.0,0.0) ,(5.0/9.0, np.sqrt(3.0/5.0)) ],
    4: [
        ( (18.0+np.sqrt(30.0))/36.0, -np.sqrt( 3.0/7.0 - 2.0/7.0*np.sqrt( 6.0/5.0 ) ) ),
        ( (18.0+np.sqrt(30.0))/36.0,  np.sqrt( 3.0/7.0 - 2.0/7.0*np.sqrt( 6.0/5.0 ) ) ),
        ( (18.0-np.sqrt(30.0))/36.0, -np.sqrt( 3.0/7.0 + 2.0/7.0*np.sqrt( 6.0/5.0 ) ) ),
        ( (18.0-np.sqrt(30.0))/36.0,  np.sqrt( 3.0/7.0 + 2.0/7.0*np.sqrt( 6.0/5.0 ) ) )
        ]
    }
def Quad_Table(i,j):
    return [ (W0GP*W1GP, Z0GP, Z1GP)
             for W0GP,Z0GP in Quad_Tables[i]
             for W1GP,Z1GP in Quad_Tables[j] ]

def GaussQuadGrid(Nrx,Nry, Nptx, Npty, start, L,H):
    start = np.array(start)

    table = Quad_Table(Nptx,Npty)
    l = float(L)/float(Nrx)
    h = float(H)/float(Nry)
    
    x = np.zeros((Nptx*Nrx*Npty*Nry,len(start)),dtype=np.double)
    w = np.zeros((Nptx*Nrx*Npty*Nry,1         ),dtype=np.double)
    for j in xrange(Nry):
        for i in xrange(Nrx):
            for k,(wz,xz,yz) in enumerate(table):
                x[j*Nrx*Npty*Nptx + Nptx*Npty*i+k,:] \
                  = [ start[0] + L*float(i)/(Nrx) + (xz+1.0)/2.0*L/(Nrx) ,
                      start[1] + H*float(j)/(Nry) + (yz+1.0)/2.0*H/(Nry) ]
                w[j*Nrx*Npty*Nptx + Nptx*Npty*i+k]  = wz * l*h/4.0
    return w,x

#
# Evaluating the RKPM method
#
def eval_rkpm(X,u, y, grad=False):
    #H,r = cf.Graphers.Build_Proximity_Graph_Given_Length(X,Ndesired, cutoff , y)
    r = SupR*np.ones(y.shape[0])
    H = cf.Graphers.Build_Proximity_Graph_Variable(X,r, y)
    #H = cf.Graphers.Build_Proximity_Graph_Uniform(X,SupR, y)
    if grad:
        ke = hr.kernel_interp_grad_u
    else:
        ke = hr.kernel_interp_u
    ui, = cf.Assemble(ke, H,
                      {'x':(X,dm_ptvec),
                       'sup':(r,dm_ptsca),
                       'u':(u,dm_ptsca),
                       'y':(y,dm_ptvec)},
                      {'uy':(dm_ptsca,)},
                      ndof=y.shape[0])
    return ui
def eval_shape(X, i, sample=100, grad=False):
    a = np.zeros(X.shape[0])
    a[i] = 1.0
    xprobe = cf.PP.init_grid(sample,sample,[-L,-L],[2*L,0],[0,2*L])
    shapei = eval_rkpm(X,a, xprobe, grad)
    return xprobe,shapei
def eval_rkpm_grid(X, ui, sample=100, grad=False):
    xprobe = cf.PP.init_grid(sample,sample,[-L,-L],[2*L,0],[0,2*L])
    u = eval_rkpm(X,ui, xprobe, grad)
    return xprobe,u

#
# Output routines
#
def write_cloud(fname, X, nodefields):
    import cornflakes as cf
    Hcloud = cf.Hypergraph()
    for l in xrange(X.shape[0]): Hcloud.Push_Edge([l])
    cf.GraphIO.write_graph(fname, Hcloud, X,nodefields)

def surfplot(Z, u, ax=None):
    import scipy.spatial
    from matplotlib import pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    tri = scipy.spatial.Delaunay(Z)
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        show=True
    else:
        show=False
    ax.plot_trisurf(Z[:,0], Z[:,1], u, triangles=tri.simplices, cmap=plt.cm.Spectral,linewidth=0, antialiased=False)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if show:
        plt.show()

def surfplotly(Z, u):
    import scipy.spatial
    import plotly as ply
    from plotly.graph_objs import Surface
    import plotly.figure_factory as FF
    ply.offline.init_notebook_mode()

    tri = scipy.spatial.Delaunay(Z)
    fig1 = FF.create_trisurf(x=Z[:,0], y=Z[:,1], z=u,
                                 #width=400,height=400,
                         simplices=tri.simplices,
                         plot_edges=False,
                         title="surface plot", aspectratio=dict(x=1, y=1, z=1))
    ply.offline.plot(fig1, filename="Torus.html")
