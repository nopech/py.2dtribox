# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 21:43:00 2018

@author: Raphi
"""
import numpy as np

class Mesh:
    nodes = []
    faces = []
    elements = []
    
    nNodes = 0
    nFaces = 0
    nElements = 0
    nBoundaryFaces = 0
    nIntFaces = 0
    
    def __init__(self, nodes, faces, elements, owners, nbFaces, nbElements):
        
        # Computations for nodes
        for n in range(len(nodes)):
            self.nodes.append(Node(nodes[n]))
        
        # Computations for faces
        for f in range(len(faces)):
            # Compute area of face (lenght of face)
            v1F = self.nodes[int(faces[f,0])].centroid
            v2F = self.nodes[int(faces[f,1])].centroid
            dF = np.subtract(v2F, v1F) # Vector subtraction
            a = np.linalg.norm(dF) # Get the magnitude of the difference
            cF = 0.5*np.add(v1F, v2F)
            # Add face with computed values to the list of faces
            if f < len(nbFaces):
                tempface = Face(faces[f], owners[f], nbFaces[f], a, cF)
            else:
                tempface = Face(faces[f], owners[f], -1, a, cF)
            self.faces.append(tempface) 
        
        # Computations for elements
        for e in range(len(elements)):
            elFaces = np.empty((0))
            elFaces = np.where(owners == e)
            elFaces = np.append(elFaces, np.where(nbFaces == e))
            # Compute the volume of a element (area of element)
            v1E = self.nodes[int(elements[e,0])].centroid
            v2E = self.nodes[int(elements[e,1])].centroid
            v3E = self.nodes[int(elements[e,2])].centroid
            dE = np.concatenate(([v1E],[v2E],[v3E])) # Concatenates all vectors
            cE = np.mean(dE, axis=0) # Compute the arithmetic mean of all vectors
            dEsq = np.pad(dE, ((0,0), (0,1)), 'constant', constant_values = (1)) # Add column of 1's at the right side so make matrix square
            vE = 0.5*np.absolute(np.linalg.det(dEsq)) # Compute the area of the element with the determinant
            self.elements.append(Element(elements[e], elFaces, nbElements[e], cE, vE))
        
        # Computation for overall mesh infos
        self.nNodes = len(nodes)
        self.nFaces = len(faces)
        self.nElements = len(elements)
        self.nBoundaryFaces = len(faces)- len(nbFaces)
        self.nIntFaces = self.nFaces - self.nBoundaryFaces

class Node:
    centroid = np.empty((0))
#    faces = np.empty((0))
#    elements = np.empty((0))
    
    def __init__(self, node):
        self.centroid = node
#        self.faces = faces
#        self.elements = elements

class Face:
    iNodes = np.empty((0))
    iOwner = 0
    iNeighbor = 0
    area = 0
    centroid = np.empty((0))
    
    def __init__(self, nodes, owner, neighbor, area, centroid):
        self.iNodes = nodes
        self.iOwner = owner
        self.iNeighbor = neighbor
        self.area = area
        self.centroid = centroid
        
class Element:
    iNodes = np.empty((0))
    iNeighbors = np.empty((0))
    iFaces = np.empty((0))
    nNeighbors = 0
    centroid = np.empty((0))
    volume = 0
    
    def __init__(self, nodes, faces, neighbors, centroid, volume):
        self.iNodes = nodes
        self.iFaces = faces
        self.iNeighbors = neighbors
        self.nNeighbors = len(self.iNeighbors)
        self.centroid = centroid
        self.volume = volume