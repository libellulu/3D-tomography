
import numpy as np
from algo_tomo_new import final_function
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

projections,maxz=final_function(30,3,0.979,0.112,30,30,3,100,8,8,8)

sensor_index=list(np.arange(0,270))
fig,axes= plt.subplots(1,len(projections[0]))
all_projections=np.zeros_like(projections[0])
for i in sensor_index:
    all_projections=np.maximum(projections[i], all_projections)
    print ('maximum of projection %d'%i,np.max(projections[i]))
for ax,cross_section in zip(axes,all_projections):
    ax.imshow(cross_section)
plt.tight_layout()
plt.show()
