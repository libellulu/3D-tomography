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


def function_creation(nb_x, nb_y , spacing,z):

    final_array=[]
    length=(nb_x-1)*spacing
    #print('length of a line of ccd:', length)


    for x in range (nb_x):

        for y in range (nb_y):
            Monarray=[-length/2+x*spacing,-length/2+y*spacing,z]
            #print('monarray',Monarray)
            final_array.append(Monarray)
    final_array=np.array(final_array)
    #print('final',final_array)
    return(final_array)

A=function_creation(9,9,2.5,0)

number_of_points_per_line=100
t=np.linspace(0,15,number_of_points_per_line)


def offset_pinhole(matrix, coordinate):
    offset_matrix=matrix+coordinate
    #print (offset_matrix)
    return offset_matrix

PDarray_position1=[0,0,106]
PDarray_position2=[0,118,0]
PDarray_position3=[0,5+115*np.cos(-np.pi*82.5/180),115*np.sin(-np.pi*82.5/180)]

new = A.copy()
new_2=A.copy()
new_forCCD3=A.copy()

new[:,2]=new_2[:,1]*-1
new[:,1]=new_2[:,2]*-1
B=new.copy()

CCD_1=offset_pinhole(A,[0,0,97])
CCD_2=offset_pinhole(B,[0,109,0])

theta= (np.pi)*173/180

#formula of a rotation on the axis x
rotation_matrix=np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])

rotate_coordinate=np.dot(rotation_matrix,new_forCCD3.T).T

CCD_3=offset_pinhole(rotate_coordinate,[0,5+102*np.cos(-np.pi*82.5/180),102*np.sin(-np.pi*82.5/180)])


fig = plt.figure()

ax = fig.add_subplot(1,2,1, projection='3d')

list_x_CCD1=[]
list_y_CCD1=[]
list_z_CCD1=[]
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
    ax.plot(x_vector_2,y_vector_2,z_vector_2)

for detector_3 in CCD_3:

    x_vector_3=(PDarray_position3[0]+t*(detector_3[0]-PDarray_position3[0]))
    y_vector_3=(PDarray_position3[1]+t*(detector_3[1]-PDarray_position3[1]))
    z_vector_3=(PDarray_position3[2]+t*(detector_3[2]-PDarray_position3[2]))

    ax.plot(x_vector_3,y_vector_3,z_vector_3,c='green')
save=(x_vector,y_vector,z_vector)

print('save',save[0][2])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

origin = np.array([0, 0, 0])
#axis and radius
p0 = np.array([-50, 0, 0])
p1 = np.array([50, 0, 0])
R = 97
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

#plt.show()

#function to have a scalar corresponding to the intensity over a line of sight (los)
def example_g(x,y,z,radius):
    if np.sqrt(x**2+y**2)>radius:
        g=0
    else:
        g=1
    return g

#this function gives several points of 3 coordinates [x,y,z], and all these points constitute
#an array that makes the Line Of sight
def LOS_creation(listx,listy,listz):
    vector_coord=[]
    for n in range(0,len(listx)):
        #print('survived')
        new_vect=[listx[n],listy[n],listz[n]]
        vector_coord.append(new_vect)
    vector_coord=np.array(vector_coord)
    return vector_coord
print('mylos',LOS_creation(list_x_CCD1[1],list_y_CCD1[1],list_z_CCD1[1])[0])

#I need to use each coordinate [x,y,z] of vector_coord inside of function_g
#to each g([x,y,z]) , it correspond a scalar value of intensity of the plasma
#then in my function integration I need to add this scalar with g[x,y,z] of the following point in my LOS_creation
# as there are 10 points [x,y,z] in each LOS the integration of ONE LOS will be a sum of 10 values
#then the function integration (or another function using it) will need to do the same sum on another LOS_creation
#example: sum on LOS_creation(list_x_CCD1[0],list_y_CCD1[0],list_z_CCD1[0]) has to be followed by :
#LOS_creation(list_x_CCD1[1],list_y_CCD1[1],list_z_CCD1[1]) and this until the final index (number of LOS in one CCD)
#then the same for the two other CCD
print('first len',len(list_x_CCD1))

#function that integrates on each point of the line, but it's a problem if the LOS is only defined by 2 points
def integration_no_interval_CCD(function_g,number_of_points_per_line,listx,listy,listz):
    intensity_list_LOS = defaultdict(list)
    remember_coord=defaultdict(list)

    for i in range(0,len(listx)):
        for n in range(0,number_of_points_per_line):
            A=LOS_creation(listx[i],listy[i],listz[i])[n]
            intensity_list_LOS[i].append(function_g(A[0],A[1],A[2],47))
            remember_coord[i].append({'x':A[0], 'y':A[1], 'z':A[2]})
    print('coordinate of the point number 60 in the line number 6 of the CCD1',remember_coord[6][60])
    print('intensity of that point', intensity_list_LOS[6][60])
    return intensity_list_LOS, remember_coord


integration_no_interval_CCD(example_g,100,list_x_CCD1,list_y_CCD1,list_z_CCD1)

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
        for n in range (0,number_of_interval):
            new_distance=new_distance+interval_size
            new_y=new_distance*np.cos(theta)*np.sin(phi)
            new_x=new_distance*np.cos(theta)*np.cos(phi)
            new_z=new_distance*np.sin(phi)
            #print('new',new_distance,new_x,new_y,new_z)
            #print('g=',function_g(new_x,new_y,new_z,47))
            intensity_list_LOS[i].append(function_g(new_x,new_y,new_z,47))
            remember_coord[i].append({'x':new_x, 'y':new_y, 'z':new_z})

    print('seconde fonction')
    print('intensity of the line 7',intensity_list_LOS[7])
    print('coord that go with the points of the line', remember_coord[7])
    return intensity_list_LOS,remember_coord



integration_with_interval(example_g,list_x_CCD1,list_y_CCD1,list_z_CCD1,50)
    # for i in range(0,len(listx)):
    #     for n in range(0,number_of_points_per_line):
    #         A=LOS_creation(listx[i],listy[i],listz[i])[n]
    #         #intensity_list_LOS[i].append(function_g(A[0],A[1],A[2],47))
    #     print('integration over line '+str(i) )
    # print('number1',intensity_list_LOS[7])
    # return intensity_list_LOS








#create empty array : arr = np.array([])
