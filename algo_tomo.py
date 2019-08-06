#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:22:10 2019

@author: lucie
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from math import sqrt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import norm
from collections import defaultdict

"""We begin here by the parameter used in the function of this program
these parameter are following the real values used in ISTTOK on the setup of
august 2019. If changes occurs, these parameter can be changed without affecting
the good use of the functions
"""

#position of the array of detectors (meaning the CCD) in the current setup
PDarray_position1=[0,0,106]
PDarray_position2=[0,118,0]
PDarray_position3=[0,5+115*np.cos(-np.pi*82.5/180),115*np.sin(-np.pi*82.5/180)]

#position of the pinholes of each cameras
Pinhole_coord1=[0,0,97]
Pinhole_coord2=[0,109,0]
Pinhole_coord3=[0,5+102*np.cos(-np.pi*82.5/180),102*np.sin(-np.pi*82.5/180)]

#create a line of point to use for the creation of the Line of sight
number_of_points_per_line=100
t=np.linspace(0,15,number_of_points_per_line)


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
def draw_cylinder(radius=100):

    """Function that will plot the cylinder representing the vessel
    This function has to be called in a plt.show()
    By default the radius of the cylinder is 100 mm, because of the ISTTTOK
    value
    """
    origin = np.array([0, 0, 0])
    #axis and radius
    p0 = np.array([-50, 0, 0])
    p1 = np.array([50, 0, 0])
    R = radius
    #vector in direction of axis
    v = p1 - p0
    #find magnitude of vector
    mag = norm(v)
    #unit vector in direction of axis
    v = v / mag
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= norm(n1)
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 100)
    theta = np.linspace(0, 2 * np.pi, 100)
    #use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)
    #generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    ax.plot_surface(X, Y, Z,color='plum')

    return X,Y,Z


#what we could call the main

#creation of the 3 CCD at the good place in space
CCD_1,CCD_2,CCD_3=creation_of_3D_sensor_in_space(function_creation(4,4,2))

#preparing the plot
fig = plt.figure()
ax = fig.add_subplot(1,2,1, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#call function for plot+show
draw_cylinder(100)
list_x_CCD3, list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1=lists_for_LOS_draw(CCD_1,CCD_2,CCD_3)
plt.show()

#function to have a scalar corresponding to the intensity over a line of sight (los)
def example_g(x,y,z,radius):
    if np.sqrt(x**2+y**2)>radius:
        g=1
    else:
        g=1
    return g

#this function gives several points of 3 coordinates [x,y,z], and all these points constitute
#an array that makes the Line Of sight
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


#I need to use each coordinate [x,y,z] of vector_coord inside of function_g
#to each g([x,y,z]) , it correspond a scalar value of intensity of the plasma
#then in my function integration I need to add this scalar with g[x,y,z] of the following point in my LOS_creation
# as there are 10 points [x,y,z] in each LOS the integration of ONE LOS will be a sum of 10 values
#then the function integration (or another function using it) will need to do the same sum on another LOS_creation
#example: sum on LOS_creation(list_x_CCD1[0],list_y_CCD1[0],list_z_CCD1[0]) has to be followed by :
#LOS_creation(list_x_CCD1[1],list_y_CCD1[1],list_z_CCD1[1]) and this until the final index (number of LOS in one CCD)
#then the same for the two other CCD


#function that integrates on each point of the line, but it's a problem if the LOS is only defined by 2 points
def integration_no_interval_CCD(function_g,listx,listy,listz):
    intensity_list_LOS = defaultdict(list)
    remember_coord=defaultdict(list)

    for i in range(0,len(listx)):
        for n in range(0,number_of_points_per_line):
            A=LOS_creation(listx[i],listy[i],listz[i])[n]
            intensity_list_LOS[i].append(function_g(A[0],A[1],A[2],47))
            remember_coord[i].append({'x':A[0], 'y':A[1], 'z':A[2]})
    print('coordinate of the point number 1 in the line number 6 of the CCD1',remember_coord[6][0])
    print('coordinate of the last point in the line number 6 of the CCD1',remember_coord[6][-1])
    print('intensity of that point 1', intensity_list_LOS[6][0])
    print('value of the integral along line number 6', np.sum(intensity_list_LOS[6]))
    return intensity_list_LOS, remember_coord


integration_no_interval_CCD(example_g,list_x_CCD1,list_y_CCD1,list_z_CCD1)

print('last elemennt',LOS_creation(list_x_CCD1[1],list_y_CCD1[1],list_z_CCD1[1])[-1])

def integration_with_interval(function_g,listx,listy,listz,interval_size):
    intensity_list_LOS = defaultdict(list)
    print('length inside',len(listx))
    remember_coord=defaultdict(list)


    for i in range(0,len(listx)):
        vector_ref_begining=LOS_creation(listx[i],listy[i],listz[i])[0]
        #print('ref',vector_ref_begining)
        vector_ref_end=LOS_creation(listx[i],listy[i],listz[i])[-1]
        #print(vector_ref_end)
        dist_x=(vector_ref_end[0]-vector_ref_begining[0])**2
        dist_y=(vector_ref_end[1]-vector_ref_begining[1])**2
        dist_z=(vector_ref_end[2]-vector_ref_begining[2])**2
        distance = np.sqrt(dist_x+dist_y+dist_z)
        #print('distance', distance)
        theta=np.arcsin((vector_ref_end[2]-vector_ref_begining[2])/distance)
        #print('theta',theta)
        #print ('le cos', np.cos(theta))
        if 0.9999<(vector_ref_end[1]-vector_ref_begining[1])/(distance*np.cos(theta))<1.001 or -1.001<(vector_ref_end[1]-vector_ref_begining[1])/(distance*np.cos(theta))<-0.999 :
            #print('here')
            phi=np.arccos((vector_ref_end[0]-vector_ref_begining[0])/(distance*np.cos(theta)))
        else:
            phi=np.arcsin((vector_ref_end[1]-vector_ref_begining[1])/(distance*np.cos(theta)))
            #print('value',(vector_ref_end[1]-vector_ref_begining[1])/(distance*np.cos(theta)))
        #print('phi',phi)
        number_of_interval=int((distance/interval_size))
        #print('numb of int', number_of_interval)
        new_distance=0
        new_y=0
        new_x=0
        new_z=106
        for n in range (0,number_of_interval):
            new_distance=new_distance+interval_size
            new_y=new_distance*np.cos(theta)*np.sin(phi)
            new_x=new_distance*np.cos(theta)*np.cos(phi)
            new_z=new_distance*np.sin(phi)
            #print('new',new_distance,new_x,new_y,new_z)
            #print('g=',function_g(new_x,new_y,new_z,47))
            intensity_list_LOS[i].append(function_g(new_x,new_y,new_z,47))
            remember_coord[i].append({'x':new_x, 'y':new_y, 'z':new_z})
    intensity_list_LOS[i].append(function_g(vector_ref_end[0],vector_ref_end[1],vector_ref_end[2],47))
#    remember_coord[i].append({'x':last, 'y':last, 'z':last})
    print('seconde fonction')
    print('coordinate of the point number 1 in the line number 6 of the CCD1',remember_coord[6][0])
    print('coordinate of the last point in the line number 6 of the CCD1',remember_coord[6][-1])
    print('intensity of that point 1', intensity_list_LOS[6][0])
    print('value of the integral along line number 6', np.sum(intensity_list_LOS[6]))
    return intensity_list_LOS,remember_coord



integration_with_interval(example_g,list_x_CCD1,list_y_CCD1,list_z_CCD1,1)
    # for i in range(0,len(listx)):
    #     for n in range(0,number_of_points_per_line):
    #         A=LOS_creation(listx[i],listy[i],listz[i])[n]
    #         #intensity_list_LOS[i].append(function_g(A[0],A[1],A[2],47))
    #     print('integration over line '+str(i) )
    # print('number1',intensity_list_LOS[7])
    # return intensity_list_LOS








#create empty array : arr = np.array([])
