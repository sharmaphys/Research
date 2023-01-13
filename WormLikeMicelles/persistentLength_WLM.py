#!/usr/bin/env python
"""
   Description:
   Calculates the persistence length of worm like micelle system
   where the end points are dynamic. First the end points of micelle
   with longest path is found. Once the two end points are found, 
   shortest path joining the two end points are found, coordinates
   of the corresponding points and bond length joining each points
   are calculated required to calculate the overall contour length
   of the micelle. Persistence length is calculated from the expo-
   nential fit of contour length vs the averages of cosines of vector
   between the first segment and rest, second and all other and so on.
  
    <cos(theta)> = exp(-S/Lp): 
   theta is the angle between the bond segments, S is the contour 
   length, and Lp is the persistence length.
                                

   Depending on the previous frame start/end points of micelle, current 
   start and end points are changed in order to keep consistent direction
   of end to end vectors, assuming that the closest frame end to end 
   vectors do not flip more than 90 degrees. 


    Example 1. Using single .gro file:
              to calculate the persitence length of a Wormlike micelle. 
              having surfactant resname OLEA and using tail bead C1, 

             ./script_name -g filename.gro -wlm "resname OLEA and name C1" -cutoff 1.2 -rg "resname OLEA"

    Example 2: Using trajectory files:

             ./script_name -f traj.xtc -s topol.tpr -n int(of surfactant) -b 19000 -e 20000 
             -wlm "resname OLEA and name C1" -cutoff 1.2 -rg "resname OLEA"

    Example 3: Calculating bond correlation function using multiple frames:

             ./script_name -f multi  -wlm "resname OLEA and name C1" -cutoff 1.2 -sep 5     

# Version 0.3                        +
"""
__author = "Hari Sharma <hari.sharma@uconn.edu>"


from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from numpy import array, sum
import numpy as np
from numpy import linalg as LA
import os
import matplotlib.pyplot as plt
from pylab import polyfit
#from IPython.display import Image
import timeit
import itertools
from itertools import chain
import structure, topology, mass_data
import dijkstra_shortest_path
from structureFactor import static_structure_factor
from structureFactor_2d import static_structure_factor_2d

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()


class Grid():
     def __init__(self, coord, box_min, box_max,  cutoff=None):

         if cutoff == None:
            self.cutoff = 0.35
         else:
            self.cutoff = cutoff
         self.coord  = coord
         # note that box_min and box_max are different than min and max of coord,
         # they are used to create the grids.
         self.box_min =  box_min
         self.box_max =  box_max
         self.box_dim = self.box_max - self.box_min 
         self.nGridsX = int(((self.box_max[0] - self.box_min[0])/self.cutoff))
         self.nGridsY = int(((self.box_max[1] - self.box_min[1])/self.cutoff))
         self.nGridsZ = int(((self.box_max[2] - self.box_min[2])/self.cutoff))

         self.ngrids = np.array([self.nGridsX, self.nGridsY, self.nGridsZ])
    
         self.invdelta = [ngrids/box for ngrids, box in zip(self.ngrids, self.box_dim)]
         #or can be used self.realcutoff as:
         #self.realcutoff = [box/ngrids for box, ngrids in zip(self.box_dim, ngrids)]

     def gen_grid(self):
         """
         generates grid with atoms assigned in grids.

         Parameters:  coord: (nx3) np array
         ----------   cutoff : float for distance cutoff, cell(grid) size
    
         returns:     grid: 3d list
         -------
         """
         # generate list in 3 dimension
         grid = [[[ [] for k in range(self.nGridsZ)] \
                       for j in range(self.nGridsY)] \
                       for i in range(self.nGridsX)]

         grid_coord = []
         for i, pos in enumerate(self.coord):
             x_pos = int((pos[0]-self.box_min[0])*self.invdelta[0])%self.nGridsX
             y_pos = int((pos[1]-self.box_min[1])*self.invdelta[1])%self.nGridsY
             z_pos = int((pos[2]-self.box_min[2])*self.invdelta[2])%self.nGridsZ
             grid[x_pos][y_pos][z_pos].append(i)
         return grid

def get_neighbor(cell, grid, nGrids):
    """
    Gives the atom indices of atoms located in neighbors of current cell.
    Each cell has 27 neighbors including itself. So, for each 
    cell, looks for atoms in current cell and atoms in neighbor cells.

    for eg: neighbors of cell [2,5,7] are :
    [1, 6, 7], [2, 6, 7], [3, 6, 7]
    [1, 5, 7], [2, 5, 7], [3, 5, 7]
    [1, 4, 7], [2, 4, 7], [3, 4, 7]
    ------------
    times 3 if you change in z dimension.    
 
    Cell is array of grid index of current atom

    Input:  cell,
    -----
    """    
    [grid_x, grid_y, grid_z] = cell
    [nGridsX, nGridsY, nGridsZ] = nGrids
    
    # index of cell at boundary
    ltEdgeX, ltEdgeY, ltEdgeZ = 0, 0, 0
    rtEdgeX, rtEdgeY, rtEdgeZ = nGridsX-1, nGridsY-1, nGridsZ-1
    
    neigh_atoms = []
    nei_cells = (-1, 0, 1)
  
    for x, y, z in itertools.product(nei_cells, nei_cells, nei_cells):
        # use of modulo implements the pbc of cells
        gx = (grid_x + x) % nGridsX   
        gy = (grid_y + y) % nGridsY
        gz = (grid_z + z) % nGridsZ
        neigh_atoms += grid[gx][gy][gz]
    return neigh_atoms


def _dist_sq(a, b):
    "a and b are array of single coordinate, size (1x3)"
    diff = a -b
    return np.dot(diff, diff) 

def _dist_sq_pbc(a, b, box):
    """
    Calculates the distance between two points applying pbc.
    a and b are array of single coordinate, size (1x3)
    box is an array with length of box in each dimension
    """
   # applying pbc for a vector between two points
    vec = b -a
    for i in range(3):
        vec[i]= vec[i] - np.rint(vec[i]/box[i])*box[i]
    return np.dot(vec, vec)

def _vec_pbc(a, b, box):
    """
    Get a vector between two points a and b applying pbc.
    a and b are array of single coordinate, size (1x3)
    box is an array with length of box in each dimension
    """
    # applying pbc for a vector between two points
    vec = b -a
    for i in range(3):
        vec[i]= vec[i] - np.rint(vec[i]/box[i])*box[i]
    return vec

def _angle_between(a, b):
    """
    Calculates the angle between two points in degrees.
    a, b are 1 d arrays
    """
    try:
        a /= np.linalg.norm(a)
        b /= np.linalg.norm(b)
    except:
        raise ValueError('Got zero vector.')
    angle = np.arccos(np.dot(a, b))*(180/np.pi)
    if np.isnan(angle):
        if (a == b).all():
            return 0.0
        else:
            return (np.pi)*(180/np.pi)
    return angle

def get_com_within_grid(coord_short_path, coord_all, gridSize):
    if not arg.gs:
       print("....where is gridcutoff size. Use argument -gs <float>")
       exit(1)
    coord_short_path = np.array(coord_short_path)
    box_min = np.min(coord_all, axis=0)
    box_max =  np.max(coord_all, axis=0)
    box_dim   = box_max-box_min

    gd = Grid(coord_all, box_min, box_max, gridSize)
    cellsize = gd.cutoff
    invdelta = gd.invdelta
    grid = gd.gen_grid()
    nGrids = gd.ngrids
   
    distcutoff_sq = gridSize*gridSize
  
    # loop over shortestPath
    com = []
    for i in range(len(coord_short_path)):
        pos = []
        neigh_atom_ndx = []
        # index of all atoms within distance of distance criteria
        cell = np.array([ int((coord_short_path[i][0]-box_min[0])*invdelta[0])%nGrids[0],  \
                          int((coord_short_path[i][1]-box_min[1])*invdelta[1])%nGrids[1],  \
                          int((coord_short_path[i][2]-box_min[2])*invdelta[2])%nGrids[2] ])
        neigh_atom_ndx = get_neighbor(cell, grid, nGrids)
        for j in range(len(neigh_atom_ndx)):
             distsq = _dist_sq_pbc(coord_short_path[i], coord_all[neigh_atom_ndx[j]], box_dim)
             # check distance criteria
             if distsq <= distcutoff_sq:
                x_j = coord_all[neigh_atom_ndx[j]][0]
                y_j = coord_all[neigh_atom_ndx[j]][1]
                z_j = coord_all[neigh_atom_ndx[j]][2]
                pos.append([x_j, y_j, z_j])
        pos = np.array(pos)
        centerOfMass = np.mean(pos, axis=0)
        com.append([centerOfMass[0],centerOfMass[1],centerOfMass[2]] )
        #        pos.append(coord_all[neigh_atom_ndx[j]])
        

#        com.append([np.mean(pos, axis= 0)])    
    return com


def mean(x):
   """"
    Returns the arithematic mean of values in 1-d array
   """
   n = len(x)
   if n < 1:
      raise ValueError('mean requires at least one data') 
   avg = sum(x)/float(n)
   return avg

def standard_dev(x):
   """
   Get the standard deviation of the 1-d array.
   """
   n = len(x)
   if n <2:
      raise ValueError(' variance requires at least two data points')
   c = mean(x)  
   variance = sum((x -c)**2)/float(n)
   sd = sqrt(variance)
   return sd 

def traj_file():
    """
    Extract .gro files from a .xtc file from initial
    time to final time and do rest of the analysis.

    """
    if not arg.index_no:
       print('\t \n What is the index number for the group you want? \n')
       exit(0)
    if not arg.topol_tpr:
       print('\t \n .tpr file needed. \n')
       exit(0)
    if not arg.begin_time:
       print('\t \n Initial time= ?? \n')
       exit(0)
    if not arg.end_time:
       print('\t \n Final time= ?? \n')
       exit(0)
 
    xtc = arg.traj_file
    tpr = arg.topol_tpr
    t1  = arg.begin_time
    t2  = arg.end_time
    n   = arg.index_no 

    command1 = "echo %s | trjconv_mpi -s %s -f %s -b %s -e %s -n -sep -pbc whole -o frame_dump_.gro > /dev/null"  %(n, tpr, xtc, t1, t2)
    subprocess.call(command1, shell=True)
    return None

def get_bond_vec(coordinates):
    """
    Get the normalized vector between two points and its norm
    given a array of coordinates.
    """
    N      = len(coordinates)
#    coordinates = np.array(coordinates)
    nBonds = N-1
    bondVector = []
    norm       = np.zeros(nBonds)
    x, y, z = 0.0, 0.0 , 0.0
    for i in range(nBonds):
         x = coordinates[i+1][0] - coordinates[i][0]
         y = coordinates[i+1][1] - coordinates[i][1]
         z = coordinates[i+1][2] - coordinates[i][2]
         norm[i] += sqrt(x**2 + y**2 + z**2)
         
         bondVector.append([x/float(norm[i]),y/float(norm[i]),z/float(norm[i])])
    
    bondVector = array(bondVector)
    return (bondVector, norm)


def _cosine_between(a, b):
    """
    Calculates the cosine between two points as vectors.
    a, b are 1 d arrays
    """
    try:
        a /= np.linalg.norm(a)
        b /= np.linalg.norm(b)
    except:
        raise ValueError('Got zero vector.')
    cosine = np.dot(a, b)
    return cosine

def persistence_ete_projection(coordinates):
    """
    Calculates the persistence length along the backbone 
    using the projection of the backbone end-to-end vector
    Re on the generic kth spring.

    lp(k) = < (rk, k+1/(|rk, k+1|)).Re >

    This is the Flory's approach and is considered to be the 
    exact and without any approximation.

    Exact persistence length is Lp = < sum over bonds ri.r1 >/ <r>, 
    where <r> is the average bond length, r1 is the first bond 
    vector and ri are all bond vectors.ie average sum of the 
    projection of backbone end to end vector on the first vector.
    """
    nBonds = len(coordinates)-1 # no of bonds
    bondVector = []
    bond_vector_normalized = []
    norm       = np.zeros(nBonds)
    x, y, z = 0.0, 0.0 , 0.0
    for i in range(nBonds):
         x = coordinates[i+1][0] - coordinates[i][0]
         y = coordinates[i+1][1] - coordinates[i][1]
         z = coordinates[i+1][2] - coordinates[i][2]
         norm[i] += sqrt(x**2 + y**2 + z**2)
         bond_vector_normalized.append([x/float(norm[i]),y/float(norm[i]),z/float(norm[i])])
         bondVector.append([x, y, z])

    coordinates = np.array(coordinates)
    # end to end vector
    re = coordinates[-1]-coordinates[0]
    # average bond length
    avg_bond_length = np.mean(norm)
    # for projection of Ree on individual bond
    lp_k= []
    lp_1 = []
    for i in range(nBonds):
        lp_k.append(np.dot(bond_vector_normalized[i], re))
        lp_1.append(np.dot(bondVector[0], bondVector[i]))

    lp_exact =  float(sum(lp_1)/avg_bond_length)

    return np.absolute(np.array(lp_k))

def vector_autocorrelation(eteVector):
    """
    Calculates the normalized end-to-end vector autocorrelation.
     <r(t).r(0)>/(r(0)**2) , where r(0) is the end to end vector at
    time t = 0 and r(t) at time t. Also returns the angular distribution. 

    input: of end to end vector, nx3 

    returns :  array of normalized autocorrelation
            :  array of angular distribution
    """
    N = len(eteVector)
    
    acf = []
    angle_acf = []
    for i in range(N):
        dotProd = np.dot(np.array(eteVector[0]), np.array(eteVector[i]))
        norm0 = np.linalg.norm(np.array(eteVector[0]))
        norm1 = np.linalg.norm(np.array(eteVector[i]))
        cos_theta = dotProd/float(norm0*norm1)   
        acf.append(cos_theta)
        theta = _angle_between(np.array(eteVector[0]), np.array(eteVector[i]))
        angle_acf.append(theta)
    return np.array(acf), np.array(angle_acf)
  

def get_segment_bond_vector(coordinates, sep):
    """
    Given the coordinates of (backbone), get the vector
    connecting between atoms separated by sep atoms.
 
    """
    N = len(coordinates)
    if sep > N:
       print("\t Atom separation greater than the number of atoms in backbone. \n")
       exit(0)
    nSegments = int(N/sep)
    bondVector = []
    for i in range(nSegments):
        ii = i * sep
        jj = (ii + sep) -1  # -1 is to get the correct indicies
        vec = coordinates[jj]-coordinates[ii]
        bondVector.append([vec[0], vec[1], vec[2]])
    return bondVector   # return as list 
    

def correlation_function_with_first(coordinates):
    """
     Calculate the Cos(theta) between the normalized vectors.
     Specifically, dot product of first vector with all the rest
     of the vectors.

     <u(0).u(i)> = <cos(theta)> 
    """
    N = len(coordinates)
    (bondVector, norm)    = get_bond_vec(coordinates)   
#    coordinates = np.array(coordinates) 
    nBonds = N-1
    contourLength = [[] for i in range(nBonds)]
    correlationFunction = [[] for i in range(nBonds)]
    
    # for correlation only with first segment.
    CL = []
    CF = []    
    s = 0
    k = 0
    cFunc = 0
    for i in range(nBonds):
        s += norm[i]
        cFunc = np.inner(bondVector[0], bondVector[i])
        CL.append(s)
        CF.append(cFunc)
    # return only positive values of CF
    CLPos = []
    CFPos = []

    for i in range(len(CL)):
        if CF[i]>0:
              CLPos.append(CL[i])
              CFPos.append(CF[i])
        else:
           break
  
#    return(CLPos, CFPos, s)
    return(CL, CF, s)


def correlation_function_all(coordinates):
    """
     Calculate the Cos(theta) between the normalized vectors.
     Specifically, dot product of first vector with all the rest
     of the vectors, dot product of second vector and all the
     rest of vectors except first vector and so on.
    """
    N = len(coordinates)
    (bondVector, norm)    = get_bond_vec(coordinates)
#    coordinates = np.array(coordinates) 
    nBonds = N-1
    contourLength = [[] for i in range(nBonds)]
    correlationFunction = [[] for i in range(nBonds)]

    for i in range(nBonds): 
        k = 0 
        s = 0                      
        cFunc = 0
        n = 1
        # to take steps of 2 for alternate bonds, as only trans bonds alligns(similar to 
        # how gromacs calculate, use step size as 2 below.
        for j in range(i, nBonds, 1):
           # get the contour length(total bond length from bond i to bond j included, (not distance from i to j))
           s = sum(norm[i:j+1])
#           cFunc = _cosine_between(np.array(coordinates[1]-coordinates[0]), np.array(coordinates[i+1]-coordinates[i]))
           cFunc += np.inner(bondVector[i], bondVector[j])
           cFunc_avg = float(cFunc/n)  
           contourLength[k].append(s)
           correlationFunction[k].append(cFunc_avg)
           k += 1
           n += 1   
    # remove emtpy lists 
    CL = [x for x in contourLength if x]
    CF = [x for x in correlationFunction if x]

    nn = len(CL)
    if len(CL)!= len(CF):
       raise ValueError('Different Length of Correlation function and Contour length')

    avgCL = [[] for i in range(nn)]
    avgCF = [[] for i in range(nn)]
    for i in range(nn):
      avgCL[i] = sum(CL[i])/float(len(CL[i]))
      avgCF[i] = sum(CF[i])/float(len(CF[i]))
 
    totalCL = sum(avgCL[-1])
    # only positive values of correlation function, assuming it should decay to 0
    avgCLPos = []
    avgCFPos = []
    
    for i in range(len(avgCL)):
        if avgCF[i]>0:
              avgCLPos.append(avgCL[i])
              avgCFPos.append(avgCF[i])
        else:
           break

#    return(avgCLPos, avgCFPos, totalCL)
    return(avgCL, avgCF, totalCL)

#

def fit_exponential_decay(x, y):
    r"""Fit a function to an exponential decay
    .. math::  y = \exp(-x/a)

    Parameters
    ----------
    x, y : array_like
      The two arrays of data
    Returns
    -------
    a : float
      The coefficient *a* for this decay
    Notes
    -----
    This function assumes that data starts at 1.0 and decays to 0.0
    Requires scipy
    """

    from scipy.optimize import curve_fit

    # get only positive values
    xPos = []
    yPos = []

    for i in range(len(x)):
        if y[i]>0:
              xPos.append(x[i])
              yPos.append(y[i])
        else:
           break

    xPos = np.array(xPos)
    yPos = np.array(yPos)
    def expfunc(x, a):
        return np.exp(-x/a)

    a = curve_fit(expfunc, xPos, yPos)[0][0]
    return a

def fit_linear(x,y):
    from scipy.optimize import curve_fit

# small matrix will simply have the [j]th position of element 1 in long matrix.
#def generate_small_matrix(matrixA, N):
#    matrixAShort = [[] for i in range(N)] # list of list with different length of row
#    for i in range(N):
#       for j in range(N):
#         if matrixA[i][j] == 1 :
#            matrixAShort[i].append(j)
#    print len(matrixAShort)
#    return matrixAShort
         

# multiply symmetric matrix with small matrix generated
#def multiply_matrix(matrixA, matrixAShort, N):
#    result = np.zeros((N, N))
#    for i in range(N):      
#        for j in range(N):
#           if matrixAShort[i]!= 0:
#             for k in range(len(matrixAShort[i])):
#                result[i][j] += matrixA[j][matrixAShort[i][k]]
#             result[j][i] = result[i][j]
#    return result 


def longest_path(coordinates, distanceCutoff):
   """
   Finds the longest path from the given point coordinates. Distance 
   cutoff is the distance within which two points are considered 
   connected.
   
   Input: n x 3 array of coordinates
          float for the cutoff of distance   
 
   Returns: int( start, end, longest step)
            adjacency matrix ( n x n array)
            number of coordnates(points)
            nxn array of actual distance between i & j th points.
   """
   
   N = len(coordinates)
   # adjacency Matrix(array)
   A = np.zeros((N, N)) 
   # distance Matrix(array)
   D = np.zeros((N, N)) 
   #store real distances between points within cutoff
   realDist = np.zeros((N,N)) 

   # D Matrix
   for i in range(N):
     for j in range(N):
       if j != i:
        D[i][j] = -1
   
   # A matrix
   for i in range(N):
     for j in range(N):
       if j > i:
         dist=np.sqrt(np.sum((coordinates[i,:] - coordinates[j,:])**2))
         if dist <= distanceCutoff:
           A[i][j] = 1
           A[j][i] = 1 
           realDist[i][j] = dist
           realDist[j][i] = dist


   currentA = A
 #  matrixAShort = generate_small_matrix(currentA, N)
   currentD = 1
   foundNewPath = True
   n = 1
   while(foundNewPath==True):
     foundNewPath = False
     for i in range(N):
       for j in range(N):
          if D[i][j] == -1 and currentA[i][j] != 0:
           foundNewPath = True
           D[i][j] = currentD
           D[j][i] = currentD
     currentD += 1
   
     if foundNewPath==True:
       n += 1
       if arg.gro_file:
         if n == 5:
             print("\nFinding the End Points of Micelle...\n")
         elif n==45:
             print("Almost there... !\n")
#       tempA = currentA
#       currentA = multiply_matrix(currentA, matrixAShort, N)
#       currentA = np.matrix(A)*np.matrix(tempA)
       currentA = LA.matrix_power(A, n)
       
## find start and end point of longest path
   longest = -1
   start = 0
   end   = 0
   for i in range(N):
      for j in range(i, N):
        if D[i][j] > longest:
           start = i
           end = j
#           if start>end:
#              start, end = end, start
           longest = D[i][j]

   return (start, end, longest, A, realDist)

        
def get_path(start, end, A, N):
    """
    Returns the shortest path that can be traversed from start to end.
    A is the adjacency matrix and N the number of points(coordinates)

    Input: int (start, end, N)
           array nx3 (N)
    returns: 1-d array of shortest path
    """
    route = [[] for i in range(N)]
    if start == end:
      route =[]
    
    else:
      visited = [ False for i in range(N)]
      visited[start] = True
      route[start].append(start)
      toCheckThisRound = []
      toCheckNextRound = []
      for i in range(N):
         if A[start][i]==1:
             FROM = start
             TO = i 
             toCheckNextRound.append([FROM, TO])
      # for each depth, search is bfs.
      while(len(toCheckNextRound) != 0):
            toCheckThisRound =  toCheckNextRound
            toCheckNextRound = []
            for list  in toCheckThisRound:
                  if visited[list[1]] != True:
                     visited[list[1]] = True
                     # for each step, keep accumulating the route, where it came from
                     # [[......] 
                     route[list[1]].append(route[list[0]])
                     route[list[1]].append(list[1])
                     for j in range(N):
                          if A[list[1]][j]==1:
                             toCheckNextRound.append([list[1],j])

     # extract connecting points of route from start to end
     # i.e, element from nested list having structure like [[[[] n1] n2] n3]
     # eg. if we go from a > b > c > d > e, final list will be [[[[a b], c], d], e]
      shortestPath = []
      while(route[end][0]!=start):
        shortestPath =  [route[end][1]] + shortestPath
        route[end] = route[end][0]  
      else:
        shortestPath =  [route[end][0]] + shortestPath
      return np.array(shortestPath)

def get_shortest_path_dijkstra(start, end, adjMat, distMat, N):
#    print('start:', start, 'end:', end)
    G = dijkstra_shortest_path.Graph()
    for i in range(N):
        G.add_vertex(str(i))
#    print("Finding Dijkstra's Shortest Path..")
    for i in range(N):
        for j in range(N):
            if j >i :
               if adjMat[i][j]==1:
                  G.add_edge(str(i), str(j), distMat[i][j])
                  G.add_edge(str(j), str(i), distMat[j][i])

    sp = dijkstra_shortest_path.shortest_path(G, str(start), str(end))
    path = []
    for i in sp:
        path.append(int(i))

    return np.array(path)

def get_com_neighbors_within(A, shortestPath, coordinates):
    """
    Calculates the com of points within cutoff distance used to 
    calculate the adjacency matrix, shortestPath.
    Input:
    ----   A = adjacency matrix nxn array
           shortestPath = 1d array 
           coordinates  = nx3 array of positions, all coordinates
    """
    n = len(shortestPath)  
    m = len(coordinates)
    com = []
    for i  in range(n):
        pos = []
        x = coordinates[shortestPath[i]][0]
        y = coordinates[shortestPath[i]][1]
        z = coordinates[shortestPath[i]][2]
        pos.append ([x, y, z])
        for j in range(m):
            if A[shortestPath[i]][j]==1:
                x_j = coordinates[j][0]
                y_j = coordinates[j][1]
                z_j = coordinates[j][2]
                pos.append([x_j, y_j, z_j])      
        pos = np.array(pos)
        centerOfMass = np.mean(pos, axis=0) 
        com.append([centerOfMass[0],centerOfMass[1],centerOfMass[2]] )
    return(com)

def get_coord_short_path(coordinates, shortestPath):
    n = len(shortestPath)
    coord_shortest = []
    for i in range(n):
        x = coordinates[shortestPath[i]][0]
        y = coordinates[shortestPath[i]][1]
        z = coordinates[shortestPath[i]][2]
        coord_shortest.append([x,y,z])

    return(coord_shortest)

def avgUnequalArray(unequalArray):
    """
    Returns nx2 array of average of columns and its standard deviation 
    from an array of arrays having unequal length.
    Example: if a = [[1],[1,2,3],[2,3,4,9]], average of array(a) returns:
             [1.33333333 2.5  3.5   9.  ]
    """
    output = []
    maxlength = 0
    for arr in unequalArray:
        if len(arr)>maxlength:
           maxlength = len(arr)
     
    for index in range(maxlength):
        temp = []
        for arr in unequalArray:
            if index < len(arr):
               temp.append(arr[index])
        if len(temp)>1:
           stdev = standard_dev(np.array(temp))
        else:
           stdev = 0 
        output.append([mean(temp), stdev])
    return np.array(output)
        
def main():
    if len(argv) < 2:
       print(parser.print_help())
       exit(1)


    if arg.gro_file:
       nFrame = 1
       grofile = arg.gro_file
       gridSize = arg.gs
 
       output = open('cf_vs_cl_wlm.xvg','w')
       output.write("@    title \"Correlation Function\" \n")
       output.write("@    xaxis  label \"Contour Length (s)(nm)\" \n")
       output.write("@    yaxis  label \" C = <u(0).u(s)> \"\n")
       output.write("@TYPE xy""\n")
    
  
       ete_correlation = open('end_to_end_correlation_wlm.xvg', 'w')
       ete_correlation.write("@    title \"Persistence Length lp(k) vs k\" \n")
       ete_correlation.write("@    xaxis  label \"k \" \n")
       ete_correlation.write("@    yaxis  label \" Persistence Length lp(k) \"\n")
       ete_correlation.write("@TYPE xy""\n")

       output3 = open('shortest_path.dat','w')
 
       results = {}
       rg = 0.0
       ete = 0.0
       cl = 0.0
       lpwlm = 0.0
       flex = 0.0


       if arg.rg:
          molecule = structure.Molecule(filename=grofile)
          selection_string1 = arg.rg
          (coord, mass) = molecule.select_atoms(selection_string1,mass=True)
          rg += topology.radius_of_gyration(coord, mw = mass)
          results['Radius of gyration'] = rg


       if arg.wlm:
          molecule = structure.Molecule(filename=grofile)
          box_size = molecule.box
          selection_string2 = arg.wlm
          coordinates = molecule.select_atoms(selection_string2)
          
          coord_all = molecule.select_atoms("resname OLEA ")
          n = len(coordinates)
          distanceCutoff = arg.cutoff
          (start, end, longest, adjMat, realDist) = longest_path(coordinates, distanceCutoff)
#          pathWalked = get_path(start, end, adjMat, n)
          pathWalked = get_shortest_path_dijkstra(start, end, adjMat, realDist, n)
#          coord_path = get_coord_short_path(coordinates, pathWalked)
          coord_path = get_com_neighbors_within(adjMat, pathWalked, coordinates)         
#          coord_path =  get_com_within_grid(pathWalked, coord_all, gridSize)
          structure.write_grofile("coord-short" ,coord_path, "AAAA", "BB", box_size)
          (CL, CF, clTotal) =  correlation_function_all(coord_path)

          cl += clTotal

          for i in range(len(CL)):
             output.write(" %8.3f    %8.3f\n" %(CL[i], CF[i]))

          lp_k = persistence_ete_projection(coord_path)
          for i in range(len(lp_k)):
              ete_correlation.write(" %8.3f    %8.3f\n" %(i+1, lp_k[i]))

          lpwlm += fit_exponential_decay(CL, CF)
          ete += topology.end_to_end_distance(np.array(coord_path[0]), np.array(coord_path[-1]))
          flex += cl/float(lpwlm)


          if arg.sf:
             sf = open('sf.xvg','w')
             sf.write("@    title \"Structure Factor\" \n")
             sf.write("@    xaxis  label \" q (1/nm)\" \n")
             sf.write("@    yaxis  label \" S(q) \"\n")
             sf.write("@    yaxis  label \" q*S(q) \"\n")
             sf.write("@TYPE xy""\n")

             q, sq = static_structure_factor(coord_path, box_size, max_q=None, binMultiplier=None)

             # write values excluding q = 0.
             for i in range( len(q)):
                sf.write(("{:.3f}" + "\t" +"{:.3f}" + "\t" + "{:.3f}" + "\n").format(q[i], sq[i], q[i]*sq[i]))

          if arg.sf2d:
             sf = open('sf2d.xvg','w')
             sf.write("@    title \"Structure Factor\" \n")
             sf.write("@    xaxis  label \" q (1/nm)\" \n")
             sf.write("@    yaxis  label \" Sx(q) \"\n")
             sf.write("@    yaxis  label \" Sy(q) \"\n")
             sf.write("@TYPE xy""\n")

             q, sxq, syq = static_structure_factor_2d(coord_path, max_q=None)

             # write values excluding q = 0.
             for i in range( len(q)):
                sf.write(("{:.3f}" + "\t" +"{:.3f}" + "\t" + "{:.3f}" + "\n").format(q[i], sxq[i], syq[i]))

          results['Contour-Length']    = cl
          results['End-to-End'] = ete
          results['Persistence Length'] = lpwlm       
          results['Flexibility'] = flex

          output3.write('Shortest Path (according to resid): '+ str(pathWalked+1) + '\n') 
          output3.write('Path Length: '+ str(len(pathWalked))) 
          print('\t Path Length :' + str(len(pathWalked)) + '\n')

       keys = list(results.keys())
       nKeys = len(list(results.keys()))
       keyLength = [(len(key)+3)*'=' for key in keys]
       output1 = open('lp-wlm.dat', 'a')
       output1.write(('    {}' + '\t\t{:>s}(nm)'*nKeys +'\n').format("Frame", *keys))
       output1.write(('    {:=^5}'+'\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       values = []
       for i in range(nKeys):
           values.append(results[keys[i]])
       output1.write(('    {:^5}'+'\t\t{:.3f}'*nKeys +'\n').format(i+1, *values))

       k = results.keys()
       v = results.values()
       merged_all = [j for i in zip(k,v) for j in i]
       print (("\n \t {:20s} : {:.3f} (nm)"*nKeys +'\n').format(*merged_all))
       print("""
       \t Note that only +ve values of averages of cosines are used to calculate persistence 
       \t length using exponential fit. \n""")

#       print ("\n \t Rough estimate of persistence Length using interpolation = %5.3f nm\n" %lp_interpolation)


    elif arg.traj_file:
       if arg.traj_file=='multi':
          pass
       else:# arg.traj_file == 'traj'
          traj_file()  
       nFrame = 0
       process = open('process.dat','w')

       results = {} 
             
 
       if arg.rg:
            results['rg']   = []
       if arg.wlm:
            results['ete']  = []          
            results['cl']   = []
            results['lp']   = []
            results['flex'] = []
            # to save from all frames
            contourLenTotal = [] 
            corrFuncTotal = []
            lpTotalAllFrame = 0

            lpk_total =[]
            lpk_std  = []
          
            # end to end vector for all frames
            eteVecTot =  []
 
            # Segmental bond vector auto correlation
            segmentVector = []

       eteVecPrevious = np.array([0, 0, 0])

       while True:
         grofile = "frame_dump_" + str(nFrame) + ".gro"
         if not path.isfile(grofile) or path.getsize(grofile) == 0:
            break


         if arg.rg:
            molecule = structure.Molecule(filename=grofile)
            selection_string1 = arg.rg
            (coord, mass) = molecule.select_atoms(selection_string1,mass=True)
            rg = topology.radius_of_gyration(coord, mw = mass)
            results['rg'].append(rg)

         if arg.wlm:
            molecule = structure.Molecule(filename=grofile)
            box_size = molecule.box
            selection_string2 = arg.wlm
            coordinates = molecule.select_atoms(selection_string2)
            n = len(coordinates)
            distanceCutoff = arg.cutoff
            (start, end, longest, adjMat, realDist) = longest_path(coordinates, distanceCutoff)
            path1 = get_shortest_path_dijkstra(start, end, adjMat, realDist, n)
            coord_path = get_com_neighbors_within(adjMat, path1, coordinates)

            # Depending on the previous frame start/end points, current start and end points are
            # changed in order to keep consistent direction of end to end vectors, assuming that the 
            # the closest frame end to end vectors do not flip more than 90 degrees. 
            
            if eteVecPrevious.all() == 0:
               eteVecPrevious = np.array(coord_path[-1]) - np.array(coord_path[0])
            eteVecCurrent = np.array(coord_path[-1]) - np.array(coord_path[0])   
            angle = _angle_between(np.array(eteVecPrevious), np.array(eteVecCurrent))
            if angle > 90:
               # flip the coordinates 
               coord_path = np.flipud(coord_path)
            #update previous ete by current ete
            eteVecPrevious = np.array(coord_path[-1]) - np.array(coord_path[0])#eteVecCurrent               

#            path1 = get_path(start, end, adjMat, n)
#            coord_path = get_com_neighbors_within(adjMat, path1, coordinates)
#            coord_path = get_coord_short_path(coordinates, path1)

            structure.write_grofile("coord-short" ,coord_path, "AAAA", "BB", box_size)

            (CL, CF, clTotal) =  correlation_function_all(coord_path)
            cl = clTotal

            contourLenTotal.append(CL)
            corrFuncTotal.append(CF)

            lpwlm = fit_exponential_decay(CL, CF)
            ete   = topology.end_to_end_distance(np.array(coord_path[0]), np.array(coord_path[-1]))
            eteVec = np.array(coord_path[-1]) - np.array(coord_path[0])
            flex  = cl/float(lpwlm)

#            print("start:",start, "end:", end, "end-to-end-vec:", eteVec)

            # for projection of end to end vector on bond
            lpk = list(persistence_ete_projection(coord_path))
            lpk_total.append(lpk)

            eteVecTot.append(eteVec)

            # for segmental autocorrelation function
            separation = arg.sep
            bondVec = get_segment_bond_vector(np.array(coord_path), separation)
            segmentVector.append(bondVec)


            results['cl'].append(cl)
            results['ete'].append(ete)
            results['lp'].append(lpwlm)
            results['flex'].append(flex)            


#         remove(grofile)
         nFrame += 1
         process.write(" Finished frame: {}, Path Length : {} \n".format(nFrame, len(path1)))
         process.flush()

       # get the average projections of end to end vector on individual segment from all frames         
       if arg.wlm:
           lpk_tot = np.array(lpk_total)
           ete_correlation = open('end_to_end_correlation_wlm.xvg', 'w')
           ete_correlation.write("@    title \"Persistence Length lp(k) vs k\" \n")
           ete_correlation.write("@    xaxis  label \"k \" \n")
           ete_correlation.write("@    yaxis  label \" Persistence Length lp(k) \"\n")
           ete_correlation.write("@TYPE xy""  \"\n")
           ete_correlation.write("@ view 0.15, 0.15, 0.75, 0.85 \"\n")
           ete_correlation.write("@ legend on  \"\n")
           ete_correlation.write("@ legend box on \"\n")
           ete_correlation.write("@ legend loctype view \"\n")
           ete_correlation.write("@ legend 0.78, 0.8  \"\n")
           ete_correlation.write("@ legend length 2  \"\n")
           ete_correlation.write("@ s0 legend \" lp(k)\" \"\n")
           ete_correlation.write("@ s1 legend \"stdev\" \n")


           lpk_avg = avgUnequalArray(lpk_tot)
           for i in range(len(lpk_avg)):
              ete_correlation.write(" %8.3f    %8.3f    %8.3f\n" %(i+1, lpk_avg[i][0], lpk_avg[i][1]))
    
         # get the average of contour length and correlation function from all frames
       if arg.wlm:
           contourLen_tot = np.array(contourLenTotal)
           corrFunc_tot = np.array(corrFuncTotal)
           cfVsCl = open('cf_vs_cl_wlm_total.xvg','w')
           cfVsCl.write("@    title \"Correlation Function\" \n")
           cfVsCl.write("@    xaxis  label \"Contour Length (s)(nm)\" \n")
           cfVsCl.write("@    yaxis  label \" C = <u(0).u(s)> \"\n")
           cfVsCl.write("@TYPE xy""\n")

           contourLen_tot_avg = avgUnequalArray(contourLen_tot)
           corrFunc_tot_avg = avgUnequalArray(corrFunc_tot)

           for i in range(len(contourLen_tot_avg)):
              cfVsCl.write(" %8.3f    %8.3f\n" %(contourLen_tot_avg[i][0], corrFunc_tot_avg[i][0]))
           lpTotalAllFrame = fit_exponential_decay(contourLen_tot_avg[:,0], corrFunc_tot_avg[:,0])


         # get the end to end vector autocorrelation function using end to end vector from all frames
       if arg.wlm:
           acf, theta  = vector_autocorrelation(eteVecTot)
           acfVsT = open('acf-wlm.xvg','w')
           acfVsT.write("@    title \"Auto Correlation Function\" \n")
           acfVsT.write("@    xaxis  label \"Time (ps)\" \n")
           acfVsT.write("@    yaxis  label \" <r(t).r(o)>/r(0)**2 \"\n")
           acfVsT.write("@TYPE xy""\n")

           for i in range(len(acf)):
              acfVsT.write(" %8.3f    %8.3f    %8.3f \n" %( i, acf[i], theta[i]))

         # get the segmental bond autocorrelation from all frames
       if arg.wlm:
            len_first_row = len(segmentVector[0])
            acf_segments_frame = np.zeros((len(segmentVector), len(segmentVector[0])))

            for i in range(len(segmentVector)):
                len_row = len(segmentVector[i])
                if len_row > len_first_row:
                   len_row = len_first_row
                for j in range(len_row):   
                    if len(segmentVector[i][j]) !=0:  # we need to check this coz the path length in wlm might 
                                                      # not be the same all the time
                       dotProd = np.dot(segmentVector[0][j], segmentVector[i][j])
                       norm0 =   np.linalg.norm(segmentVector[0][j])
                       norm1 =   np.linalg.norm(segmentVector[i][j])
                       cos_theta = dotProd/float(norm0*norm1)
                       p2 = 0.5*(3*(cos_theta)**2 -1)
                       acf_segments_frame[i][j] +=p2

            acf_segment_avg = np.mean(acf_segments_frame, axis = 1)

            seg_acf = open('segacf-wlm.xvg', 'w')
            seg_acf.write("@    title \"Segmental auto correlation\" \n")
            seg_acf.write("@    xaxis  label \"time(t) \" \n")
            seg_acf.write("@    yaxis  label \" P2 \"\n")
            seg_acf.write("@TYPE xy""\n")
#            for i in range(len(acf_segment_avg)):
#                seg_acf.write(" %5d    %8.3f \n" %(i, acf_segment_avg[i]))
#            seg_acf.close()
            # write the first column as the mean of all the other columns            
            for i in range(len(acf_segments_frame)):
#                for j in range(len(acf_segments_frame[0])):
                    seg_acf.write((" {:>3d}" + "{:>8.3f}" +   "{}" +"\n").format(i, acf_segment_avg[i], " ".join(format(x, "8.3f") for x in acf_segments_frame[i])))
            seg_acf.close()

       keys = list(results.keys())
       nKeys = len(keys)
       keyLength = [(len(key)+3)*'=' for key in keys]
       output1 = open('lp-wlm.dat', 'w')
       output1.write(('    {}' + '\t\t{:>s}(nm)'*nKeys +'\n').format("Frame", *keys))
       output1.write(('    {:=^5}' + '\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       for i in range(len(results[keys[0]])):
           values = []
           for j in range(nKeys):
               values.append(results[keys[j]][i])
           output1.write(("    {:^5d}" + "\t\t{:.3f}"*nKeys + "\n").format(i+1, *values))

     
       # get the standard deviation of data   
       mean_of_keys = []
       std_of_keys= []
       for i in range(nKeys):
           avgKeys = mean(np.array(results[keys[i]]))
           stdKeys = standard_dev(np.array(results[keys[i]]))
           mean_of_keys.append(avgKeys)
           std_of_keys.append(stdKeys)
 
       # get merged list with average and sd in sequence
       merged_values = [j for i in zip(mean_of_keys, std_of_keys) for j in i] 
       output1.write(('    {:=^5}' + '\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       output1.write(('    {}' + '\t {:.3f} +/- {:.3f}'*nKeys +'\n').format("Avgerages: ", *merged_values))

       # get merged keys with corresponding values
       merged_all = []
       for i in range(len(keys)):
           merged_all.append(keys[i])
           for j in range(i*2, (i*2)+2):
               merged_all.append(merged_values[j])

       for i in range(len(merged_all)):
           if merged_all[i]=='rg':
              merged_all[i] = "Radius of Gyration"
           elif merged_all[i] == "ete":
              merged_all[i] = "End-to-End distance"
           elif merged_all[i] == "lp":
              merged_all[i] = "Persistence Length"
           elif merged_all[i] == "cl":
              merged_all[i] = "Contour Length"
           elif merged_all[i] == "flex":
              merged_all[i] = "Flexibility"

       print (("\n \t {:20s} : {:.3f} +/- {:.3f}"*nKeys +'\n').format(*merged_all))
       print("""
       \t Note that only +ve values of averages of cosines are used to calculate persistence 
       \t length using exponential fit. \n
       \t Total Persistence Length = {:.3f}""".format(lpTotalAllFrame))
       
       output1.write("\t Total Persistence Length from all frame = {:.3f}".format(lpTotalAllFrame))


    end_time = timeit.default_timer()
    totalTimeTaken = end_time-start_time
    stdout.write("\t Time Taken = %i Seconds.\n \n" %totalTimeTaken)

def args():
    parser = ArgumentParser(description=__doc__, 
			   formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-f","--traj_file", help="""pass "multi" to use multiple .gro files for analysis.
                         or pass .xtc(trajectory) file for analysis.""",type=str)

    parser.add_argument('-g','--gro_file', help=""" Single .gro file for analysis. """, type=str)

    parser.add_argument('-b','--begin_time', default= 1,  help=""" Starting time for analysis.\n """, type=int)
   
    parser.add_argument('-e','--end_time', help="""  End time for analysis. """, type=int)

    parser.add_argument('-n','--index_no', help=""" Index number for polymer chain. """, type=int)

    parser.add_argument('-s', '--topol_tpr', default="topol.tpr", help=""" topol.tpr file. """, type=str)

    parser.add_argument('-rg', help = """ 
                        Calculate Rg of selected atoms (eg. resname ST1 ST2 and name C31 C21 C11) .""", type= str)

    parser.add_argument('-wlm', help = """ 
                        Calculate persistence length of coordinates that gives the longest path between 
                        end points of worm like micellar system. Atom selection example (eg. resname OLEA and name O2).""", type= str)


    parser.add_argument('-cutoff', help = """ 
                        Distance Cutoff to calculate distance and adjacency matrix.""", type= float)

    parser.add_argument('-gs', help = """ 
                        Dimension of grid size  .""", type= float)

    parser.add_argument('-sf', help = """ 
                        Calculates the static structure factor of the selected coordinates. """, action='store_true')

    parser.add_argument('-sf2d', help = """ 
                        Calculates the 2d static structure factor. """, action='store_true')
 
    parser.add_argument('-sep', dest= 'sep' , default = '3', help = """ 
                        Separation of atoms to calculate the segmental bond vector autocorrelation.
                        Default value of separation is 3. 
                        eg. -sep 5 calculates bond between 1st and 5th, 6th and 10th atom and so on.
                        Writes segacf.xvg file with second column the mean of segment correlation and 
                        all the rest with individual segment correlation.
                        """, type=int)


    return parser, parser.parse_args() 

if __name__ =='__main__':
   parser, arg = args()
   main()
