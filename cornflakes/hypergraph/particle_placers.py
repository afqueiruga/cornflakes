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
                          + random.uniform(-dev,dev,len(start))
    return x
