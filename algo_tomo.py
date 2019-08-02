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
print('PD3', PDarray_position3)


new = A.copy()
new_2=A.copy()
new_forCCD3=A.copy()

new[:,2]=new_2[:,1]*-1
new[:,1]=new_2[:,2]*-1
B=new.copy()
#print( 'A:',A)
#print('B',B)


CCD_1=offset_pinhole(A,[0,0,97])
CCD_2=offset_pinhole(B,[0,109,0])
#print(CCD_2)

theta= (np.pi)*173/180

#formula of a rotation on the axis x
rotation_matrix=np.array([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])

rotate_coordinate=np.dot(rotation_matrix,new_forCCD3.T).T
print('rotate coorrd', rotate_coordinate)
CCD_3=offset_pinhole(rotate_coordinate,[0,5+102*np.cos(-np.pi*82.5/180),102*np.sin(-np.pi*82.5/180)])
print('CCD3', CCD_3)

fig = plt.figure()

ax = fig.add_subplot(1,2,1, projection='3d')

####for n in range(1,2):


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
    #print('x', x_vector_3)
    y_vector_3=(PDarray_position3[1]+t*(detector_3[1]-PDarray_position3[1]))
    #print('y',y_vector_3)
    z_vector_3=(PDarray_position3[2]+t*(detector_3[2]-PDarray_position3[2]))
    #print('z',z_vector_3)
    ax.plot(x_vector_3,y_vector_3,z_vector_3,c='green')






ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');




z_c=0
y_c=0
R=47

x_centers=np.linspace(-35,35,30)
thetas=np.linspace(0, 6.28,20)
for x in x_centers:
    y_list=np.cos(thetas)*R-y_c
    z_list=np.sin(thetas)*R+z_c
    x_list = np.ones(20)*x
    ax.plot(x_list,y_list,z_list, c='plum',alpha=0.7)




plt.show()



#create empty array : arr = np.array([])
