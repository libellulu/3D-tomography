import numpy as np

theta= np.pi /2

rotation_matrix=np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
# mat=np.zeros((3,3)) to eeclare 3 by 3 offset_matrix
#mat= np.zeros() to declare an empty final_array
#mat=np.zeros(3) one array of 3 columns

print(rotation_matrix)

vector=np.array([[0,2,1],[0,2,1],[0,2,1],[3,5,0]])
#to transpose a vector : vector.T
rotate_coordinate=np.dot(rotation_matrix,vector.T).T
print(rotate_coordinate)
#np.dot => normal matrix product
