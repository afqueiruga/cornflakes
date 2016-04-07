import numpy as np
import numpy.random as random

def init_grid(Nx,Ny, start,e1,e2, dev=0.0):
    start = np.array(start)
    e1=np.array(e1)
    e2=np.array(e2)
    x = np.zeros((Nx*Ny,len(start)),dtype=np.double)
    for j in xrange(Ny):
        for i in xrange(Nx):
            x[j*Nx+i,:] = start + float(i/(Nx-1.0))*e1 + float(j/(Ny-1.0))*e2 \
                           + (0.0 
                             if i==0 or i==Nx-1 or j==0 or j==Ny-1else
                          random.uniform(-dev,dev,len(start)))
    return x
def init_cube(Nx,Ny,Nz, start,e1,e2,e3, dev=0.0):
    start = np.array(start)
    e1=np.array(e1)
    e2=np.array(e2)
    e3=np.array(e3)
    x = np.zeros((Nx*Ny*Nz,len(start)),dtype=np.double)
    for k in xrange(Nz):
        for j in xrange(Ny):
            for i in xrange(Nx):
                x[k*Nx*Ny+j*Nx+i,:] = start + float(i/(Nx-1.0))*e1 \
                              + float(j/(Ny-1.0))*e2 \
                              + float(k/(Nz-1.0))*e3 \
                          + (0.0 
                             if i==0 or i==Nx-1 or j==0 or j==Ny-1 or k==0 or k==Nz-1 else
                          random.uniform(-dev,dev,len(start)))
    return x

def sphere_test(cen,rad):
    return lambda x: ((x-cen).dot(x-cen)<rad)

def select_nodes(X,fil):
    return np.asarray(np.where(map(fil,X))[0],dtype=np.intc)


def delete_particles(X, conditions):
    marked = np.ones(X.shape[0],dtype=np.bool)
    for i,a in enumerate(X):
        for c in conditions:
            if c(a):
                marked[i] = False
                break
    return X[np.where(marked)]
    
