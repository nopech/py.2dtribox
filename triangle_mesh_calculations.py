# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 20:18:36 2018

@author: Raphi
"""

import numpy as np
import matplotlib.pyplot as plt
import Mesh
from scipy.spatial import Delaunay

###############################################################################
# Mesh parameters
plot_mesh = True
L = 1
Nx = 2
Ny = 2
hx = L/Nx
hy = L/Ny
points = np.zeros(((Nx+1)*(Ny+1),2))

###############################################################################
# Generate points
i = 0
for x in np.arange(0, L+hx, hx):
    for y in np.arange(0, L+hy, hy):
        points[i][0] = x
        points[i][1] = y
        i = i+1

###############################################################################
# Generate mesh  
tri = Delaunay(points)
plt.triplot(points[:,0], points[:,1], tri.simplices)
plt.plot(points[:,0], points[:,1], 'o')

if plot_mesh:
    for j, p in enumerate(points):
        plt.text(p[0]-0.003, p[1]+0.003, j, ha='right') # label the points
    for j, s in enumerate(tri.simplices):
        p = points[s].mean(axis=0)
        plt.text(p[0], p[1], '#%d' % j, ha='center') # label triangles
    plt.show()

###############################################################################
# Create list of faces with index connectivity
faces = np.empty((0,2))
owners = np.empty((0))
neighbors = np.empty((0))
boundaries = np.empty((0,2))
boundaryowners = np.empty((0))

# Loop through all simplices
for simplex in range(tri.nsimplex):
    
    # Generate faces arount current simplex
    simp_faces = np.empty((0, 2))
    for vertex in range(len(tri.simplices[simplex])):
        if vertex < len(tri.simplices[simplex])-1:
            simp_faces = np.concatenate((simp_faces, [[tri.simplices[simplex, vertex], tri.simplices[simplex, vertex+1]]]))
        else:
            simp_faces = np.concatenate((simp_faces, [[tri.simplices[simplex, vertex], tri.simplices[simplex, 0]]]))
    
    # Loop through all faces of current simplex
    for face in range(len(simp_faces)):
        
        # Shift index because the neighbor is on the opposite side of the index
        if face+2 < len(simp_faces)-1:
            i = face +2
        else:
            i = face-1 
            
        # If neighbor num is greater than simplex num, add face to list of faces and put simplex num in the list of owners
        # and add neighbor num to the list of neighbors.
        # If neighbor num is -1, add simplex num to boundaryowners and add face to the list of boundary faces
        if tri.neighbors[simplex, i] > simplex:  
            faces = np.concatenate((faces, [simp_faces[face]]))
            owners = np.append(owners, simplex)
            neighbors = np.append(neighbors, tri.neighbors[simplex, i])
        elif tri.neighbors[simplex, i] < 0:
            boundaries = np.concatenate((boundaries, [simp_faces[face]]))
            boundaryowners = np.append(boundaryowners, simplex)

# Append boundary faces to the list of faces
nboundaries = len(boundaries)
faces = np.concatenate((faces, boundaries))
owners = np.append(owners, boundaryowners)

# Put mesh data in a useful datastructure           
m = Mesh.Mesh(nodes = points,  
                   faces = faces, 
                   elements = tri.simplices, 
                   owners = owners, 
                   nbFaces = neighbors, 
                   nbElements = tri.neighbors
                   )

###############################################################################
# blubb