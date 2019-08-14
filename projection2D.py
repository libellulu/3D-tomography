
from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import MultiLineString, LineString
from mpl_toolkits import mplot3d
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import norm
from collections import defaultdict


def function_creation(nb_x, nb_y , spacing,z=0):
    """This function returns an array that represents a 3D matrix for a new
    tomography sensor.

    Parameters
    ----------
    nb_x : integer
    number of line in x to make a x*y grid for the sensor.
    The same goes for nb_y
    spacing : integer
    The space between the centers of two cell, the size of the cell
    z: integer
    originally height distance between pinhole and PD array. Has been fixed to
    zero because of further calculation and offset specific to each CCD

    Return
    ------
    final_array : Array of arrays
    Gives the final matrix created for one CCD. It is an array of (x,y,z)
    coordinates. This coordinates are the ones of the cells of the matrix.
    """

    final_array=[]
    length=(nb_x-1)*spacing

    for x in range (nb_x):

        for y in range (nb_y):
            Monarray=[-length/2+x*spacing,-length/2+y*spacing,z]
            final_array.append(Monarray)
    final_array=np.array(final_array)

    return(final_array)
def offset_pinhole(matrix, coordinate):
    """Function that calculate the new place of the detector in the space
    The detector is a matricx of coordinate (x,y,z) that are offset through
    this function

    Parameters
    ----------
    matrix : an array of array
    It is the output of function creation , an array of (x,y,z) coordinate
    coordinates : array
    the coordinates [x,y,z] of the point that we offset from.
    Here we offset from the pinhole coordinates

    Return
    ------
    offset_matrix : an array of array
    an array of the new coordinate
    """
    offset_matrix=matrix+coordinate
    return offset_matrix
def lists_for_LOS_draw(CCD_1,CCD_2,CCD_3):
    """Function that we use to draw the cones of tomography.
    It calculates list of x points, of y , and of z, that are used to draw
    the lines of sight, and will also remember them in a list.
    This function has to be called in a plt.show(), or followed by it
    to see the cones

    Parameters
    ----------
    CCD_1 (same for 2 and 3): array of array
    the coordinates calculated , then offset, and replaced by the
    function creation_of_3D_sensor_in_space

    Return
    ------
    list_x_CCD1: array of list
    successive lists of x coordinates that are in each LOS of CCD1
    will be the same for y, z and for CCD2 and 3
    """

    list_x_CCD1=[]
    list_y_CCD1=[]
    list_z_CCD1=[]
    list_x_CCD2=[]
    list_y_CCD2=[]
    list_z_CCD2=[]
    list_x_CCD3=[]
    list_y_CCD3=[]
    list_z_CCD3=[]
    for detector in CCD_1:
        #print('detectors:', detector)
        x_vector=-(PDarray_position1[0]+t*(detector[0]-PDarray_position1[0]))
        y_vector=-(PDarray_position1[1]+t*(detector[1]-PDarray_position1[1]))
        z_vector=(PDarray_position1[2]+t*(detector[2]-PDarray_position1[2]))
        list_y_CCD1.append(y_vector)
        list_x_CCD1.append(x_vector)
        list_z_CCD1.append(z_vector)
        ax.plot(x_vector,y_vector,z_vector,c='red')

    for detector_2 in CCD_2:
        #print('detectors2:', detector_2)
        x_vector_2=-(PDarray_position2[0]+t*(detector_2[0]-PDarray_position2[0]))
        y_vector_2=(PDarray_position2[1]+t*(detector_2[1]-PDarray_position2[1]))
        z_vector_2=-(PDarray_position2[2]+t*(detector_2[2]-PDarray_position2[2]))
        list_y_CCD2.append(y_vector_2)
        list_x_CCD2.append(x_vector_2)
        list_z_CCD2.append(z_vector_2)
        ax.plot(x_vector_2,y_vector_2,z_vector_2)

    for detector_3 in CCD_3:
        x_vector_3=(PDarray_position3[0]+t*(detector_3[0]-PDarray_position3[0]))
        y_vector_3=(PDarray_position3[1]+t*(detector_3[1]-PDarray_position3[1]))
        z_vector_3=(PDarray_position3[2]+t*(detector_3[2]-PDarray_position3[2]))
        list_y_CCD3.append(y_vector_3)
        list_x_CCD3.append(x_vector_3)
        list_z_CCD3.append(z_vector_3)

        ax.plot(x_vector_3,y_vector_3,z_vector_3,c='green')
    return list_x_CCD3,list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1

def LOS_creation(listx,listy,listz):
    """This function gives several points of 3 coordinates [x,y,z],
    and all these points constitute an array that makes the Line Of sight

    Parameters
    ----------
    listx: list of arrays
        Each array contains the x coordinates of each point in a given Line of
        sight in  one single CCD.
        Same goes for listy and listz.

    Returns
    -------
    vector_coord : ndarray
        This list has size [LOS* npoints *3]. LOS is the number of lines of
        sight. npoint , the number of point per LOS. 3 corresponds to
        coordinate [x,y,z]
    """
    vector_coord=[]
    for n in range(0,len(listx)):
        #print('survived')
        new_vect=[listx[n],listy[n],listz[n]]
        vector_coord.append(new_vect)
    vector_coord=np.array(vector_coord)
    return vector_coord
def creation_of_3D_sensor_in_space(matrix_already_done):
    """ Function creating my sensor, the place of those are relative to the
    center 0,0,0 of the vessel
    Requires to already have used function_creation previously to have an
    original matrix of points

    Parameters
    ----------
    matrix_already_done : the output of function_creation

    Returns
    ------
    CCD_1(same for CCD_2, CCD_3) : Array of arrays
    Gives the final matrix created for one CCD, in its new space.
    It is an array of (x,y,z) coordinates.
    These coordinates are the ones of the cells of the matrix.
    """
    #creation of 3 copies of the original matrix to use after
    new = matrix_already_done.copy()
    new_2=matrix_already_done.copy()
    new_forCCD3=matrix_already_done.copy()

    # multiply the coordinates by -1 and exchange value ofaxes y and z to rotate
    # the second CCD
    new[:,2]=new_2[:,1]*-1
    new[:,1]=new_2[:,2]*-1
    second_matrix=new.copy()

    CCD_1=offset_pinhole(matrix_already_done,Pinhole_coord1)
    CCD_2=offset_pinhole(second_matrix,Pinhole_coord2)

    theta= (np.pi)*173/180

    #Computation with formula of a rotation on the axis x of the place and angle
    #of the third camera (called CCD3) in the current setup

    #general calcul to make a rotation along x
    rotation_matrix=np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])
    rotate_coordinate=np.dot(rotation_matrix,new_forCCD3.T).T

    CCD_3=offset_pinhole(rotate_coordinate,Pinhole_coord3)

    return CCD_1,CCD_2,CCD_3

#position of the array of detectors (meaning the CCD) in the current setup
PDarray_position1=[0,0,106]
PDarray_position2=[0,118,0]
PDarray_position3=[0,5+115*np.cos(-np.pi*82.5/180),115*np.sin(-np.pi*82.5/180)]

#position of the pinholes of each cameras
Pinhole_coord1=[0,0,97]
Pinhole_coord2=[0,109,0]
Pinhole_coord3=[0,5+102*np.cos(-np.pi*82.5/180),102*np.sin(-np.pi*82.5/180)]

#preparing the plot
fig = plt.figure()
ax = fig.add_subplot(1,2,1, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')


#create a line of point to use for the creation of the Line of sight
number_of_points_per_line=100
t=np.linspace(0,15,number_of_points_per_line)

n_rows = 15  # y-axis pixel resolution
n_cols = 15  # x-axis pixel resolution

x_min = -100.
x_max = +100.


y_min = -100.
y_max = +100.


def transform(x, y):
    j = int((x-x_min)/(x_max-x_min)*n_cols)
    i = int((y_max-y)/(y_max-y_min)*n_rows)
    return (i, j)

# -------------------------------------------------------------------------

x_grid = np.linspace(x_min, x_max, num=n_cols+1)
y_grid = np.linspace(y_min, y_max, num=n_rows+1)

grid = []

for x in x_grid:
    grid.append([(x, y_min), (x, y_max)])

for y in y_grid:
    grid.append([(x_min, y), (x_max, y)])

grid = MultiLineString(grid)

# -------------------------------------------------------------------------
CCD_1,CCD_2,CCD_3=creation_of_3D_sensor_in_space(function_creation(4,4,2))
list_x_CCD3, list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1=lists_for_LOS_draw(CCD_1,CCD_2,CCD_3)

projections = []

for n in range (0,len(list_x_CCD3)):
    A=LOS_creation(list_x_CCD3[n], list_y_CCD3[n],list_z_CCD3[n])
    line = LineString([(A[0][0], A[0][1]), (A[-1][0], A[-1][1])])
    print('linestring', line)
    projection = np.zeros((n_rows, n_cols))
    for segment in line.difference(grid):
        xx, yy = segment.xy
        x_mean = np.mean(xx)
        y_mean = np.mean(yy)
        (i, j) = transform(x_mean, y_mean)
        projection[i,j] = segment.length
    projections.append(projection)
    #projections.append(projection * row.etendue)

projections = np.array(projections)

print('projections:', projections.shape, projections.dtype)

# -------------------------------------------------------------------------

fname = 'projections.npy'
print('Writing:', fname)
np.save(fname, projections)

# -------------------------------------------------------------------------

fig, axes=plt.subplots(4,4)
for i in range(16):
    axes.flatten()[i].imshow(projections[i])
plt.show()

# fig, ax = plt.subplots(ni, nj, figsize=figsize)
# for i in range(ni):
#     for j in range(nj):
#         k = i*nj + j
#         ax[i,j].imshow(projections[k], vmin=None, vmax=None)
#         ax[i,j].set_axis_off()
#
# fig.suptitle('projections (vertical camera)')
# plt.show()
#
# # -------------------------------------------------------------------------
#
# vmin = 0.
# vmax = np.max([projections[k] for k in df[df['camera']=='horizontal'].index])
#
# fig, ax = plt.subplots(ni, nj, figsize=figsize)
# for i in range(ni):
#     for j in range(nj):
#         k = i*nj + j + ni*nj
#         ax[i,j].imshow(projections[k], vmin=vmin, vmax=vmax)
#         ax[i,j].set_axis_off()
#
# fig.suptitle('projections (horizontal camera)')
# plt.show()
