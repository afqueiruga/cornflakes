from cornflakes import *
import numpy as np
from matplotlib import pylab as plt

Nsidex, Nsidey = 50,50
W = 10.0
H = 10.0

X = np.vstack([
    PP.init_grid(Nsidex,Nsidey,[-W/2.0, -H/2.0],[W,0.0],[0.0,H],0.0),
    PP.init_grid(Nsidex,Nsidey,[-W/2.0, -H/2.0],[0.5*W,0.0],[0.0,0.5*H],0.0)
    ])

Y = cflib.Remove_Duplicate_Particles_np(X, 0.1, 1.1*W/Nsidex)

from cornflakes import Graphers
H,r = Graphers.Build_Proximity_Graph_Given_Length(Y,8, 2.1*W/Nsidex)

from circles import circles
circles(Y[:,0],Y[:,1],r, fc='none')
plt.plot(Y[:,0],Y[:,1],'o')
plt.show()
