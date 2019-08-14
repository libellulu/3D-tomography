import geometry as geo
import numpy as np
line=geo.Line(0,0,0,1,1,1)
voxel=geo.Voxel(0,0,0,2,2,2)
intersection= geo.intersect(voxel,line)
print(intersection.length)
geo.Line( )
geo.Voxel()
