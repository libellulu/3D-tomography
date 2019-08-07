import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import sqrt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(1,2,1, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

x, y = np.mgrid[-1.0:1.0:30j, -1.0:1.0:30j]
# Need an (N, 2) array of (x, y) pairs.
xy = np.column_stack([x.flat, y.flat])

mu = np.array([0.0, 0.0])

sigma = np.array([.025, .025])
covariance = np.diag(sigma**2)

z = multivariate_normal.pdf(xy, mean=mu, cov=covariance)

# Reshape back to a (30, 30) grid.
z = z.reshape(x.shape)

#plt.plot(x,y)
#plt.plot(y,z)
#plt.show()

def synthetic_plasma_profile(sig,xcenter,ycenter,zcenter):
    N=1/((2*np.pi)**(3/2)*sig**3)
    x=np.linspace(0,10,10)
    y=np.linspace(0,10,10)
    z=np.linspace(0,10,10)
    g_of_xyz=N*np.exp(((x-xcenter)**2+(y-ycenter)**2+(z-zcenter)**2)/(2*sig**2))
