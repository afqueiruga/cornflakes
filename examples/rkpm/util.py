import numpy as np
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


def surfplot(Z, u):
    import scipy.spatial
    from matplotlib import pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    tri = scipy.spatial.Delaunay(Z)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.plot_trisurf(Z[:,0], Z[:,1], u, triangles=tri.simplices, cmap=plt.cm.Spectral,linewidth=0, antialiased=False)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()
def write_cloud(fname, X, nodefields):
    import cornflakes as cf
    Hcloud = cf.Hypergraph()
    for l in xrange(X.shape[0]): Hcloud.Push_Edge([l])
    cf.GraphIO.write_graph(fname, Hcloud, X,nodefields)

def surfplotly(Z, u):
    import scipy.spatial
    import plotly as ply
    from plotly.graph_objs import Surface
    import plotly.figure_factory as FF
    ply.offline.init_notebook_mode()

    tri = scipy.spatial.Delaunay(Z)
    fig1 = FF.create_trisurf(x=Z[:,0], y=Z[:,1], z=u,
                         simplices=tri.simplices,
                         title="Torus", aspectratio=dict(x=1, y=1, z=0.3))
    ply.offline.iplot(fig1, filename="Torus")
