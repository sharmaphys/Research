#! /usr/bin/env python3.6

"""
   Calculates the curvature of curved lipid surface,
   mean, gaussian and principal curvature, surface normals,
   vmd visualization arrows, 
"""

__author = "Hari Sharma"

import numpy as np
import math
import scipy.spatial as spatial
import structure, topology
from sys import argv,stdout
import timeit
from scipy import spatial
from collections import OrderedDict
#from sklearn.neighbors import NearestNeighbors
#numpy.set_printoptions(threshold=numpy.nan)
start_time = timeit.default_timer()


def draw_arrow(v1, v2, r, color):
    """ 
       Draws arrow from v1 to v2 for vmd visualization.
       v1 and v2 are arrays of shape (1,3)
    """
    v1 = v1*10 # conversion to angstrom
    v2 = v2*10
    v = (v2 -v1)*0.8 + v1 # arbitary number 0.8
    arrow = ""
    arrow += "draw color " + color + " \n" 
    arrow += "draw cylinder  " + "{" + str(v1[0]) +" "+ str(v1[1]) + " "+ str( v1[2]) + "}" + \
     "  " + "{" + str(v[0]) +" "+ str(v[1]) + " "+ str( v[2]) + "}" + " radius " + str(r) + "\n"

    arrow += "draw color " + color + " \n"
    arrow += "draw cone  " + "{" + str(v[0]) +" "+ str(v[1]) + " "+ str( v[2]) + "}" + \
     "  " + "{" + str(v2[0]) +" "+ str(v2[1]) + " "+ str( v2[2]) + "}" + " radius " + str(r*3) + "\n"
    return arrow 
      

def write_vmd_arrows(fname,head_xyz, tail_xyz):
    """
       Write the arrows to a .tcl file to load in vmd
  
    Input:
    ------
    fname : filename
    head_xyz: coordinates of head beads, np array of shape (nx3)
    tail_xyz: coordinates of head beads, np array of shape (nx3)

    output: file (filename)
    -------
    """

    f = open(fname,'w')
    for i in range(len(head_xyz)):
        f.write(draw_arrow(head_xyz[i], tail_xyz[i], 0.2, "red"))
    f.close()
    

def normal_test(head_xyz, tail_xyz):
    """
    """
    normal = []
    for i in range(len(head_xyz)):
        x = head_xyz[i][0]-tail_xyz[i][0]
        y = head_xyz[i][1]-tail_xyz[i][1]
        z = head_xyz[i][2]-tail_xyz[i][2]
        normal.append([x,y,z])
    return normal

def nearest_neighbor(coord, radius):
   """
      Finds all the points within certain radius of points.
      Uses scipy library.
   
      Inputs: coord: np array of shape (nx3)
      ------  radius: float, distance within which to search
      
      Returns: nd array of nearest neighbours  
      ------
 
      NOTE::  that sanity check has carried out that is required for
      ----    surface fit. ie neighbour list has to be at least 5 so
              that the point itself makes 6 data points required later
              for quadric fit. If not every iteration is carried out
              by increasing the radius by 0.1 until 5 neighbours are 
              found.
 
              Better to implement later to find exactly k nearest 
              neighbour rather than neighbour within certain radius
              for this purpose.
   """

   n = len(coord)
   coord_tree = spatial.cKDTree(coord)
   
   neigh_found = []
   for i in range(n):
       neigh_found.append(coord_tree.data[coord_tree.query_ball_point(coord[i,:], radius)])
   
   neigh_len = [len(neigh_found[i]) for i in range(n)]
 
   while(any(neigh_len[i]<5 for i in range(len(neigh_len)))):
        neigh_len[:] = []
        neigh_found[:] = []
        result = "\t |** Current Radius for NN Searching: %4.2f. Found insufficient to find 5 NN. |" %radius
        print( result)
        radius += 0.1
        print("\t |** Using %4.2f Radius for NN searching. "%radius + " "*(len(result)-43)+"|")
        print("\t"+" " +"-"*(len(result)-2))
        for i in range(n):
            neigh_found.append(coord_tree.data[coord_tree.query_ball_point(coord[i,:], radius)]) 
        
        neigh_len = [len(neigh_found[i]) for i in range(n)]
   return np.array(neigh_found)


def k_nearest_neighbor(coord):
   """
    Finds k nearest neighbour to a point.
   """
   knn= NearestNeighbors(n_neighbors=6)
   knn.fit(coord)
   NearestNeighbors(algorithm='auto', leaf_size=30, n_neighbors=5, p=2, \
         radius=1.0, warn_on_equidistant=True)  #p= 2 is euclidian distance

   #perform queries, to get neighbors of coord[0], first data point:
   a = knn.kneighbors(coord[0], return_distance=False)
   return a

#   pass



def surface_patch(coord, found_neigh):
    """
    Creates surface patch array associated with each points 
    and its nearest neighbour. First coordinate(row) will be the centered
    point itself and all the rest are its neighbours within certain radius.

    Input: coord: (nx3) np array 
    ------ found_neigh: nearest neighbours, array of arrays

    Returns: array of lists of dimension (n1x3), (n2x3)....
             where n1 comprises (n1-1) nearest neighbour of 1st coordinate(row)
    """
    n = len(coord)
    patch = [[] for i in range(n)]
    for i in range(n):
        patch[i].append([coord[i][0], coord[i][1], coord[i][2]])
        for row in found_neigh[i]:#range(len(found_neigh[i])):
            patch[i].append([row[0], row[1], row[2]])    


#    orderedPatch = []
#    for i in range(len(patch)):
#     oPatch = map(list, OrderedDict.fromkeys(map(tuple, patch[i])).keys())
#     orderedPatch.append(oPatch)
    return np.array(patch)
    
def fit_quadric(coord):
   """
      Fits the quadric surface:

              z = Ax**2 + By**2 + Cxy + Dx + Ey + F

      using the least square fit method. i. e. the goal is to find the 
      coefficients A, B, C, D, E, and F that gives the best fit surface
 
      In linear algebra, the problem is to solve for x,

                            Ax = b,                      
               
               [ x1**2   y1**2   x1*y1  x1  y1  1  ]       [A]       [z1]
               | x2**2   y2**2   x2*y1  x2  y2  1  |       |B|       |z2]
               | ...     ....    ....   ..  ..  .  |       |C|    =  |..]  
               |				   |       |D|       |..]
               |                                   |       |E|       |..]
               [ xm**2   ym**2   xm*ym  xm  ym  1  ](mxn)  [F](nx1)  [..] (mx1)

      where A is a matrix (mxn) m > n, (n is the number of coeffs)
      x is a vector (nx1) and b is a vector (mx1), that minimizes the norm(residuals)
      || AX -b||.
      Note that the implementation shown above is how it is used in numpy.
      
      However, more general theory is, if you need to do least square fit of 
      for example y = mx + c, for a given set of samples {(x_i, y_i)}, i = 1 to n,
      we need to determine m and c so that the line y = mx+c best fits the sample
      in the sense that the squared errors between y_i and the line values mx_i+c
      is minimized. Error is only in y direction.
      
      so, if error E(m,c) = sum(i to n)[(mx_i+c)-y_i]**2. This function if non
      negative and its graph is a paraboloid whose vertex occurs when the gradient
      satisfies del(E) = (0,0).
      
         ie. (0,0) = del(E) = 2*sum(i = 1 to n)[(mx_i +c)-y_i](x_i, 1)
      
      ==>>[ sum(i = 1 to n)(x_i**2)   sum(i = 1 to n)(x_i)] [m]  = [ sum(i = 1 to n)(x_i*y_i)]
          [ sum(i = 1 to n)(x_i)      sum(i = 1 to n)(1)]   [c]  = [ sum(i = 1 to n)(y_i)]
          
        Solving this matrix problem provides the least squares solution y = mx+c.


      Input: ndarray (mx3) array 
      -----

      Returns: 1-d array of coefficients
      -----

   """  
   x = coord[:,0]
   y = coord[:,1]
   z = coord[:,2]

   n = len(x)
   x_sq = np.power(x,2)   
   y_sq = np.power(y,2)
   xy  = x*y
   
   c1 = np.reshape(x_sq, (n, 1))
   c2 = np.reshape(y_sq, (n, 1))
   c3 = np.reshape(xy, (n, 1))
   c4 = np.reshape(x, (n, 1))
   c5 = np.reshape(y, (n, 1))
   c6 = np.ones((n, 1))
  
   b = np.reshape(z,(n,1))
   A = np.hstack((c1, c2, c3, c4, c5, c6))
   coeffs = np.linalg.lstsq(A,b)[0]

   return coeffs

def transform(coord, basis):
   """
   Transforms the coordinate from one frame to another.
   For example if vector u = (u1, u2, u3) is in standard cartesian
   coordinate system with bases (e1, e2, e3), e1 = (1, 0, 0)
   e2 = (0 , 1, 0) and e3 = (0, 0, 1). ie. u = u1*e1 + u2*e2 + u3*e3:

   then the same vector u in different frame of reference say u' with 
   bases say (p1, p2, p3) can be obtained as::
                 
                    u' = inverse([M])*u
   Matrix M is the change of basis matrix given as [ column(p1), column(p2), column(p3)]

  eg, B = {u, v}, B' = {u',v'}. Suppose the basis vectors u' and v' have following
  coordinates relative to the basis B: [u']B = [a , b].T  [v']B = [c, d].T

                           ==>  u' = au + bv
                                v' = cu + dv
  the change of coordinats matrix from B' to B is P = [[a, b], [c,d]].T
                              [V]B = P[V]B'
  ie if we know the coordinates of v relative to the basis B', then this vector in basis
  B can be obtained by simply multiplying this vector by change of coordinate matrix.
   
  
   Parameters: coord: to transform, np array of shape (nx3)
   ------      basis: to which to transform; np array of shape (3x3)
                      for 3d transformation

   Returns: transformed coordinate; np array of shape (nx3) 
   ------
                                                                                             
   """
   trans_coord = []   
   invMat = np.linalg.inv(basis)
   # matrix multiplication using numpy
   for i in range(len(coord)):
      trans_coord.append(invMat.dot(coord[i]))
#      trans_coord.append([np.matmul(invMat,coord[i])])
      
   return np.array(trans_coord)

def back_transform(coord, basis):
   """
   Transform the coordinates in Principal axes frame back to box frame

   Parameters: coord: to transform, np array of shape (nx3)
   ------      basis: to which to transform; np array of shape (3x3)
                      for 3d transformation

   Returns: transformed coordinate; np array of shape (nx3) 
   ------
                                                                                  |   |   |  
   [u]box = P[u']principal frame, P is the basis matrix with bases on columns ie [p1, p2, p3]
			                                                          |   |   |
   """

   back_trans_coord = []
   for i in range(len(coord)):
       back_trans_coord.append(np.matmul(basis,coord[i]))

   return np.array(back_trans_coord)

def smooth(coord, fitPara):
    """
     Get smoothed coordinates.
    """
    smoothed = []
    n = len(coord)
    
    for i in range(n):
        x = coord[i][0]
        y = coord[i][1]
        z = fitPara[0]*x**2 + fitPara[1]*y**2 + fitPara[2]*x*y + fitPara[3]*x + fitPara[4]*y + fitPara[5] 
        smoothed.append([x,y,z])
    return np.array(smoothed)
   

def main():
    
    molecule = structure.Molecule(filename= argv[1])
    coord_head    = molecule.get_coord_atom("NC3")
    coord_tail    = molecule.get_coord_atom("C4A","DPPC")
    neighbor = nearest_neighbor(coord_head, 1.4)

    patch     = surface_patch(coord_head, neighbor)

    # shift each patch to patch with central point on origin
    structure.write_grofile("out", patch[0], "DPPC", "NC3")

    """
    patchShifted = np.array(patch[0])-np.array(patch[0])[0]
    paxis  = topology.principal_axes(patchShifted)
    transformedPatch = transform(patchShifted, paxis.T)


    

    fit_para = fit_quadric(transformedPatch)
    smoothed_coord = smooth(transformedPatch, fit_para)
    backTr = back_transform(smoothed_coord, paxis.T) + np.array(patch[0])[0]
    print(backTr)

    structure.write_grofile("smoothed", backTr, "DPPC", "NC3")

#    backTr2 = back_transform(transformedPatch, paxis.T)

#    fit_para_direct = fit_quadric(backTr2)
#    smoothed_coord_direct = smooth(backTr2, fit_para_direct)

#    structure.write_grofile("direct", smoothed_coord_direct, "DPPC", "NC3" )


    """
    smoothed_head = []
    for i in range(len(patch)):
         x = np.array(patch[i])[0]
         patch[i] = np.array(patch[i]) -np.array(patch[i])[0]
#         structure.write_grofile("out", patch[i], "DPPC", "NC3" )

         paxis  = topology.principal_axes(patch[i])
         transformedPatch = transform(patch[i], paxis.T)
         fit_para = fit_quadric(transformedPatch)
         smoothed_coord = smooth(transformedPatch, fit_para)
         backTr = back_transform(smoothed_coord, paxis.T) 
         
         smoothed_head.append(backTr[0]+x)
    smoothed = np.array(smoothed_head)

    structure.write_grofile("smoothed_direct", smoothed, "DPPC", "NC3" )

   
    ## iterate n times for smoothing
    
    for j in range(4):
         neighbor1 = nearest_neighbor(smoothed, 1.4)
         patch1     = surface_patch(smoothed, neighbor1)
         smoothed = []
         for i in range(len(patch1)):

            x1 = np.array(patch1[i])[0]
            patch1[i] = np.array(patch1[i]) -np.array(patch1[i])[0]
            paxis1  = topology.principal_axes(patch1[i])
            transformedPatch1 = transform(patch1[i], paxis1.T)
            fit_para1 = fit_quadric(transformedPatch1)
            smoothed_coord1 = smooth(transformedPatch1, fit_para1)
            backTr1 = back_transform(smoothed_coord1, paxis1.T)

            smoothed.append(backTr1[0]+x1)
         smoothed = np.array(smoothed)
         if j ==3:
             structure.write_grofile("smoothed_head-multiple", smoothed, "DPPC", "NC3" )
    
      

#    normal = normal_test(coord_head, coord_tail)

#    structure.write_grofile("head", coord_head, "DPPC", "NC3")
#    write_vmd_arrows("arrow.tcl", coord_tail,coord_head)
    
    return None

if __name__== "__main__":
   main()                                                              
