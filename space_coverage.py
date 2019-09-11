
import numpy as np
from algo_tomo_new import final_function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

#file that shows the space covered bythe selected lines of sight.
#This example shows 180 LOS

projections,maxz=final_function(20,3,0.65,0.112,50,50,3,100)

sensor_index=list(np.arange(0,180))
fig,axes= plt.subplots(1,len(projections[0]))
all_projections=np.zeros_like(projections[0])
for i in sensor_index:
    all_projections=np.maximum(projections[i], all_projections)

for ax,cross_section in zip(axes,all_projections):
    ax.pcolormesh(np.linspace(-100,100,len(cross_section[0])+1),np.linspace(100,-100,len(cross_section)+1),cross_section)

    ax.set_aspect('equal', adjustable='box')
    ax.add_patch(mpl.patches.Circle((0, 0), 85, color='black', fill=False,linewidth=3))

plt.show()
