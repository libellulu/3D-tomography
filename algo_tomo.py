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
import geometry as geo

"""We begin here by the parameter used in the function of this program
these parameter are following the real values used in ISTTOK on the setup of
august 2019. If changes occurs, these parameter can be changed without affecting
the good use of the functions
"""

#position of the array of detectors (meaning the CCD) in the current setup
#actually the distance is pinhole to the array and not the reverse
PDarray_to_pinhole1=[0,9,0]
PDarray_to_pinhole2=[9,0,0]
PDarray_to_pinhole3=[13*np.cos(-np.pi*82.5/180),13*np.sin(-np.pi*82.5/180),0]

#position of the pinholes of each cameras
Pinhole_coord1=[5,97,0]
Pinhole_coord2=[109,0,0]
Pinhole_coord3=[5+102*np.cos(-np.pi*82.5/180), 102*np.sin(-np.pi*82.5/180),0]

#create a line of point to use for the creation of the Line of sight
number_of_points_per_line=10
t=np.linspace(0,25,number_of_points_per_line)

radius_tokamak=100


def function_creation(nb_x, nb_z , spacing_x, spacing_z, y=0):
    """This function returns an array that represents a 3D matrix for a new
    tomography sensor.

    Parameters
    ----------
    nb_x : integer
    number of line in x to make a x*z grid for the sensor.
    The same goes for nb_z
    spacing : integer
    The space between the centers of two cell, the size of the cell
    y: integer
    originally height distance between pinhole and PD array. Has been fixed to
    zero because of further calculation and offset specific to each CCD

    Return
    ------
    final_array : Array of arrays
    Gives the final matrix created for one CCD. It is an array of (x,y,z)
    coordinates. This coordinates are the ones of the cells of the matrix.
    """

    final_array = []
    length_x = (nb_x-1) * spacing_x
    length_z = (nb_z-1) * spacing_z

    for x in range(nb_x):
        for z in range(nb_z):
            Monarray=[-length_x/2+x*spacing_x, y, -length_z/2 + z*spacing_z]
            final_array.append(Monarray)
    final_array=np.array(final_array)
    print(final_array)

    return(final_array)


def offset_pinhole_and_array(matrix, coordinate_pinhole,distance_to_array):
    """Function that calculate the new place of the detector in the space
    The detector is a matrix of coordinate (x,y,z) that are offset through
    this function

    Parameters
    ----------
    matrix : an array of array
    It is the output of function creation , an array of (x,y,z) coordinate
    coordinates_pinhole : array
    the coordinates [x,y,z] of the point that we offset from.
    Here we offset from the pinhole coordinates, and we will add to this the
    distance between the pinhole and the pd array
    distance_to_array: coordinates [x,y,z], contributes to the offset

    Return
    ------
    offset_matrix : an array of array
    an array of the new coordinate
    """
    offset_matrix=matrix+coordinate_pinhole+distance_to_array
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

    #multiply the coordinates by -1 and exchange value of axes y and x to rotate
    #the second CCD
    new[:,1]=new_2[:,0]*-1
    new[:,0]=new_2[:,1]*-1
    second_matrix=new.copy()
    # theta2=np.pi/2
    # rotation_matrix_x=np.array([[1,0,0],[0,np.cos(theta2),np.sin(theta2)*(-1)],[0,np.sin(theta2),np.cos(theta2)]])
    # second_matrix=np.dot(rotation_matrix_x,new_2.T).T

    CCD_1=offset_pinhole_and_array(matrix_already_done,Pinhole_coord1,PDarray_to_pinhole1)
    CCD_2=offset_pinhole_and_array(second_matrix,Pinhole_coord2,PDarray_to_pinhole2)

    theta= -(np.pi)*172.5/180

    #Computation with formula of a rotation on the axis x of the place and angle
    #of the third camera (called CCD3) in the current setup

    #general calcul to make a rotation along x
    rotation_matrix=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    rotate_coordinate=np.dot(rotation_matrix,new_forCCD3.T).T

    CCD_3=offset_pinhole_and_array(rotate_coordinate,Pinhole_coord3,PDarray_to_pinhole3)
    return CCD_1,CCD_2,CCD_3


def lists_for_LOS_draw(CCD_1, CCD_2, CCD_3, plot_list=[]):
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
        x_vector=(detector[0]+t*(Pinhole_coord1[0]-detector[0]))
        y_vector=(detector[1]+t*(Pinhole_coord1[1]-detector[1]))
        z_vector=(detector[2]+t*(Pinhole_coord1[2]-detector[2]))
        list_y_CCD1.append(y_vector)
        list_x_CCD1.append(x_vector)
        list_z_CCD1.append(z_vector)
        if 'CCD1' in plot_list:
            ax.plot(x_vector,y_vector,z_vector,c='red')


    for detector_2 in CCD_2:
        #print('detectors2:', detector_2)
        x_vector_2=(detector_2[0]+t*(Pinhole_coord2[0]-detector_2[0]))
        y_vector_2=(detector_2[1]+t*(Pinhole_coord2[1]-detector_2[1]))
        z_vector_2=(detector_2[2]+t*(Pinhole_coord2[2]-detector_2[2]))
        list_y_CCD2.append(y_vector_2)
        list_x_CCD2.append(x_vector_2)
        list_z_CCD2.append(z_vector_2)
        if 'CCD2' in plot_list:
            ax.plot(x_vector_2,y_vector_2,z_vector_2)

    for detector_3 in CCD_3:
        x_vector_3=(detector_3[0]+t*(Pinhole_coord3[0]-detector_3[0]))
        y_vector_3=(detector_3[1]+t*(Pinhole_coord3[1]-detector_3[1]))
        z_vector_3=(detector_3[2]+t*(Pinhole_coord3[2]-detector_3[2]))
        list_x_CCD3.append(x_vector_3)
        list_y_CCD3.append(y_vector_3)
        list_z_CCD3.append(z_vector_3)
        if 'CCD3' in plot_list:
            ax.plot(x_vector_3,y_vector_3,z_vector_3, color='blue')

    plt.ylim((-600, 600))
    plt.xlim((-600, 600))
    ax.set_zlim((-600, 600))

    print(np.linalg.norm([list_x_CCD3[0][-1]-list_x_CCD3[0][0],list_y_CCD3[0][-1]-list_y_CCD3[0][0]]))
    print(np.linalg.norm([list_x_CCD3[-1][-1]-list_x_CCD3[-1][0],list_y_CCD3[-1][-1]-list_y_CCD3[-1][0]]))

    return list_x_CCD3,list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1


def find_furthest_z(listx,listy,listz):
    """Function that finds the furthest point in z among the LOS of one CCDs.

    Parameters
    ----------
    listx : list of integer
    the list of x coordinatesof all the points constituting all the LOS
    of one CCD. The return of the lists_for_LOS_draw functions
    listy, listz : idem in the other axes

    Returns
    -------
    maxx : integer
    the furthest point in the space in the z direction
    """
    list_new_z_a=[]

    for i in range(0,len(listx)):
        # print('nouveau')
        # print('i=', i)
        # print('ma liste', listx[i])
        # print('length of list i', len(listx[i]))
        # print('list i de 2', listx[i][2])
        for n in range (0, len(listx[i])):
            new_z_a=listz[i][n][np.sqrt(listx[i][n]**2+listy[i][n]**2)<100]
            if new_z_a.size >0:
                list_new_z_a.append(new_z_a)
    maxx=np.max(list_new_z_a)

    return maxx
def max_z_among_all_CCD(listx1,listy1,listz1,listx2,listy2,listz2,listx3,listy3,listz3):
    """Function that returns the furthest point in the z axes among all the LOS
     of all the CCD.

     Parameters
     ----------
     listx1 : list of integer
     list of the x coordinates of all the points constituting the
     several LOS of CCD1.
     listy1, listz1 : idem in the other axes
     listx2,listx,3: idem for the two other CCD

     returns
     -------
     supermax : integer
     the furthest z point among all CCD
    """
    max_for_CCD1=find_furthest_z(listx1,listy1,listz1)
    max_for_CCD2=find_furthest_z(listx2,listy2,listz2)
    max_for_CCD3=find_furthest_z(listx3,listy3,listz3)

    if max_for_CCD1>max_for_CCD2:
        supermax=max_for_CCD1
    else:
        supermax=max_for_CCD2
        if max_for_CCD3>supermax:
            supermax=max_for_CCD3

    return supermax
def voxel_creation(max_found,nb_voxel_x, nb_voxel_y,nb_voxel_z,radius_tokamak):
    """Function that discretize the space in voxel . The space in the one where
    LOS are present inside of the vessel.

    Parameters
    ----------
    max_found : integer
    the coordinates z of the furthest point in the z axes among all
    LOS in all CCD. Return of the function max_z_among_all_CCD.
    nb_voxel_x :integer
    number of voxel in the x axes
    same for nb_voxel_y and z
    radius_tokamak:integer
    the radius specified in the beginning of the program

    returns
    -------
    voxellist: list of the type geo.voxel. Can display only begining and end
    with "name of the voxel".start or .end
    list of voxel discretising the space

    """
    voxellist=[]
    discr_of_x= np.linspace(-radius_tokamak,radius_tokamak,nb_voxel_x+1)
    discr_of_y= np.linspace(radius_tokamak,-radius_tokamak,nb_voxel_y+1)
    discr_of_z= np.linspace(-max_found,max_found,nb_voxel_z+1)
    z0=discr_of_z[0]
    x0=discr_of_x[0]
    y0=discr_of_y[0]
    for z in discr_of_z[1:]:
        for y in discr_of_y[1:]:
            for x in discr_of_x[1:]:
                voxellist.append(geo.Voxel(x0,y0,z0,x,y,z))
                x0=x.copy()
            x0=discr_of_x[0]
            y0=y.copy()
        y0=discr_of_y[0]
        z0=z.copy()

    return voxellist
def intersection_all_LOS(matrix_voxel_created,listx,listy,listz):

    """
    """
    #first, adaptation of the LOS
    remember_line=[]
    remember_intersection=[]
    for i in range(0,len(listx)):
        point1=LOS_creation(listx[i],listy[i],listz[i])[0]
        #print('premier',point1)
        point2=LOS_creation(listx[i],listy[i],listz[i])[-1]
        #print('second',point2)
        line=geo.Line(point1[0],point1[1],point1[2],point2[0],point2[1],point2[2])
        remember_line.append(line)
        #now the intersection
        for n in range(0,len(matrix_voxel_created)):
            current_voxel=matrix_voxel_created[n]
            intersection= geo.intersect(current_voxel,line)
            remember_intersection.append(intersection.length)
            #print('intersection of the LOS'+str(i)+'with the voxel number'+str(n), intersection.length)

    return remember_intersection

def draw_cylinder(radius_tokamak):

    """Function that will plot the cylinder representing the vessel
    This function has to be called in a plt.show()
    By default the radius of the cylinder is 100 mm, because of the ISTTTOK
    value
    """
    origin = np.array([0, 0, 0])
    #axis and radius
    p0 = np.array([-20, 0, 0])
    p1 = np.array([20, 0, 0])
    R = radius_tokamak
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
    Z, X, Y = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    ax.plot_surface(X, Y, Z,color='plum')

    return X,Y,Z
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
    example : LOS_creation(list_x_CCD1[5], list_y_CCD1[5],list_z_CCD1[5])will
    return a array of all the point and their coordinate [x,y,z] constituting
    the LOS 5 of the CCD1
    """
    vector_coord=[]
    for n in range(0,len(listx)):
        #print('survived')
        new_vect=[listx[n],listy[n],listz[n]]
        vector_coord.append(new_vect)
    vector_coord=np.array(vector_coord)
    return vector_coord
#function to have a scalar corresponding to the intensity over a line of sight (los)
def example_g(x,y,z,radius=47):
    if np.sqrt(x**2+y**2)>radius:
        g=1
    else:
        g=1
    return g
def integration_LOS(function_g,listx,listy,listz):
    """Function that integrates the intensity of each point of a given line
    of sight

    Parameters
    ----------
    function_g: the fucntion of density of plasma. It itself takes 3 Parameters
    that are coordinate x,y,z and this function return the intensity in that
    point
    listx: a list of all the x coordinates for all the points
    listy, listz: same as listx
    These lists are used for LOS_creation inside of the function

    Returns
    -------
    intensity_list_LOS: a list of integers value of the intensity in each points
    of coordinates (x,y,z)
    example : intensity_list_LOS[A][B] with A being the number of the LOS in
    the CCD we chose and [B] being the number of the point in this line. This
    will return the intensity of point number B

    remember_coord: a list of coordinates for each point of each LOS
    example : remember_coord[A][B] (see above) will return the (x,y,z)
    coordinate of point number B in LOS number A
    """
    intensity_list_LOS = defaultdict(list)
    remember_coord=defaultdict(list)

    for i in range(0,len(listx)):
        vector_ref_begining=LOS_creation(listx[i],listy[i],listz[i])[0]
        vector_ref_end=LOS_creation(listx[i],listy[i],listz[i])[-1]
        dist_x=(vector_ref_end[0]-vector_ref_begining[0])**2
        dist_y=(vector_ref_end[1]-vector_ref_begining[1])**2
        dist_z=(vector_ref_end[2]-vector_ref_begining[2])**2
        length_of_line = np.sqrt(dist_x+dist_y+dist_z)
        interval_dl= length_of_line/number_of_points_per_line

        for n in range(0,number_of_points_per_line):
            A=LOS_creation(listx[i],listy[i],listz[i])[n]
            intensity_list_LOS[i].append(function_g(A[0],A[1],A[2])*interval_dl)
            remember_coord[i].append({'x':A[0], 'y':A[1], 'z':A[2]})
    return intensity_list_LOS, remember_coord
def example_of_use_of_integration():
    """Function that shows what can the function integration_LOS do
    Has to be called in a print
    """
    intensity_list_LOS, remember_coord=integration_LOS(example_g,list_x_CCD1,list_y_CCD1,list_z_CCD1)
    print('coordinate of the point number 1 in the line number 6 of the CCD1',remember_coord[6][0])
    print('coordinate of the last point in the line number 6 of the CCD1',remember_coord[6][-1])
    print('intensity of that point 1', intensity_list_LOS[6][0])
    print('value of the integral along line number 6', np.sum(intensity_list_LOS[6]))#what we could call the main
#function not working NOW , on its way to work one day
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

    print('seconde fonction')
    print('coordinate of the point number 1 in the line number 6 of the CCD1',remember_coord[6][0])
    print('coordinate of the last point in the line number 6 of the CCD1',remember_coord[6][-1])
    print('intensity of that point 1', intensity_list_LOS[6][0])
    print('value of the integral along line number 6', np.sum(intensity_list_LOS[6]))
    return intensity_list_LOS,remember_coord

#creation of the 3 CCD at the good place in space
nb_cell_x=10
nb_cell_z=3
CCD_1,CCD_2,CCD_3=creation_of_3D_sensor_in_space(function_creation(nb_cell_x, nb_cell_z, spacing_x=2, spacing_z=2))

#preparing the plot
fig = plt.figure()
ax = fig.add_subplot(1,2,1, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#call function for plot+show
draw_cylinder(100)
list_x_CCD3, list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1=lists_for_LOS_draw(CCD_1,CCD_2,CCD_3,plot_list=['CCD1','CCD3'])
#plt.show()
A=max_z_among_all_CCD(list_x_CCD3, list_y_CCD3,list_z_CCD3,list_x_CCD2,list_y_CCD2,list_z_CCD2,list_x_CCD1,list_y_CCD1,list_z_CCD1)
print("maximum z", A)

nb_voxel_x = 20
nb_voxel_y = 20
nb_voxel_z = 3

S=voxel_creation(A,nb_voxel_x,nb_voxel_y,nb_voxel_z,radius_tokamak)

list_of_intersection=[]
list_of_intersection.append(intersection_all_LOS(S,list_x_CCD1,list_y_CCD1,list_z_CCD1))
list_of_intersection.append(intersection_all_LOS(S,list_x_CCD2,list_y_CCD2,list_z_CCD2))
list_of_intersection.append(intersection_all_LOS(S,list_x_CCD3,list_y_CCD3,list_z_CCD3))
projections=np.array(list_of_intersection).flatten().reshape((nb_cell_x * nb_cell_z * 3, nb_voxel_z, nb_voxel_y, nb_voxel_x))
np.save('projections.npy',projections)

# for projection_cube in projections:
#     fig, axes = plt.subplots(1, len(projection_cube))
#     for ax, p in zip(axes, projection_cube):
#         ax.imshow(p)
