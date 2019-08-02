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


def function_creation(nb_x, nb_y , spacing,z):

    final_array=[]
    length=(nb_x-1)*spacing
    print('length of a line of ccd:', length)


    for x in range (nb_x):

        for y in range (nb_y):
            Monarray=[-length/2+x*spacing,-length/2+y*spacing,z]
            #print(Monarray)
            final_array.append(Monarray)
    final_array=np.array(final_array)
    #print(final_array)
    return(final_array)

A=function_creation(9,9,2.5,0)


t=np.linspace(0,15,10)


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

for detector in CCD_1:
    #print('detectors:', detector)
    x_vector=-(PDarray_position1[0]+t*(detector[0]-PDarray_position1[0]))
    y_vector=-(PDarray_position1[1]+t*(detector[1]-PDarray_position1[1]))
    z_vector=(PDarray_position1[2]+t*(detector[2]-PDarray_position1[2]))
    ax.plot(x_vector,y_vector,z_vector,c='red')

for detector_2 in CCD_2:
    #print('detectors2:', detector_2)
    x_vector_2=-(PDarray_position2[0]+t*(detector_2[0]-PDarray_position2[0]))
    y_vector_2=(PDarray_position2[1]+t*(detector_2[1]-PDarray_position2[1]))
    z_vector_2=-(PDarray_position2[2]+t*(detector_2[2]-PDarray_position2[2]))
    ax.plot(x_vector_2,y_vector_2,z_vector_2)

for detector_3 in CCD_3:
    #print('detectors2:', detector_2)
    x_vector_3=(PDarray_position3[0]+t*(detector_3[0]-PDarray_position3[0]))
    y_vector_3=(PDarray_position3[1]+t*(detector_3[1]-PDarray_position3[1]))
    z_vector_3=(PDarray_position3[2]+t*(detector_3[2]-PDarray_position3[2]))
    ax.plot(x_vector_3,y_vector_3,z_vector_3,c='green')


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');

origin = np.array([0, 0, 0])
#axis and radius
p0 = np.array([-50, 0, 0])
p1 = np.array([50, 0, 0])
R = 48.5
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
ax.plot_surface(X, Y, Z)

plt.show()



#create empty array : arr = np.array([])
