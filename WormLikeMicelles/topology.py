#! /usr/bin/env python3
"""
   Various analyses.
"""

__author = "Hari Sharma"

import numpy as np
from numpy import linalg as LA
import math
import structure
from sys import argv,stdout
import timeit
from scipy import spatial
#numpy.set_printoptions(threshold=numpy.nan)
start_time = timeit.default_timer()


def center_of_mass(coord):
    """
    Calculates the center of mass of the input coordinates.
    Assumes equal mass for all atoms.
 
    Parameters
    ------------

    coordinates: ndarray, shape = (natoms, 3)

    Returns
    -----------
    Center of Mass: ndarray, shape = (1, 3) ie. [x,y,z] 
    
    """
    com = np.mean(coord, axis = 0)
    return com

def moment_of_inertia(coord):
   """
   Computes the moment of inertia of the given structure(coordinates)
   relative to the center of mass .Assumes equal mass for all atoms. 

        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2));
        # Iyy = sum(m_i*(x_i^2+z_i^2));
        # Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)

   Parameters
   -------------
   Input: coord: nd array, shape (nx3)

   Returns
   -------------
   moment of inertia tensor: 3x3 array
   """
   
   com = np.mean(coord, axis = 0)
   mi = 1.0 
  #translate the coordinates to body center of mass

   pos = coord-com 

   inertiaTensor = np.zeros((3, 3), dtype=np.float64)
   # xx
   inertiaTensor[0][0] = (mi*(pos[:, 1] ** 2 + pos[:, 2] ** 2)).sum()
   # xy & yx
   inertiaTensor[0][1] = inertiaTensor[1][0] = -(mi*(pos[:, 0] * pos[:, 1])).sum()
   # xz & zx
   inertiaTensor[0][2] = inertiaTensor[2][0] = -(mi* (pos[:, 0] * pos[:, 2])).sum()
   # yy
   inertiaTensor[1][1] = (mi*(pos[:, 0] ** 2 + pos[:, 2] ** 2)).sum()
   # yz + zy
   inertiaTensor[1][2] = inertiaTensor[2][1] = -(mi*(pos[:, 1] * pos[:, 2])).sum()
   # zz
   inertiaTensor[2][2] = (mi*(pos[:, 0] ** 2 + pos[:, 1] ** 2)).sum()

   
   #inertiaTensor = np.dot(pos.transpose(), pos)
   
   return inertiaTensor    
   
def radius_of_gyration(coord, **kwargs):
   """
   Computes the radius of gyration of the input structure. 
   So far, mass is assumed to be equal 
   
   Parameters
   -------------
   Input: pos(coordinates of atoms)

   kwargs:
         mw = mass weighted, if mass weighted, you need mass 
         data(array ) of all atoms in selected group.
   
   Returns
   -------------
   Radius of gyration
      
   """
   com = np.mean(coord, axis = 0)

   # mass weighted
   mw = 'mw' in kwargs.keys()
   if mw:
      m = kwargs['mw']
   else: m = [1]*len(coord)
   total_mass = sum(m)

   #translate the com to origin   
   recenteredCom = coord -com
   N = len(recenteredCom)  #number of atoms

   # (sum(m_i*(R_i-R)^2))/sum(m_i)
   # for all atom of equal mass, Rg_sq =sum(R_i-R)^2/N, R = center 

   rgSquared = float(np.sum(m*np.sum(recenteredCom**2, axis=1))/total_mass)
   return np.sqrt(rgSquared)

def end_to_end_distance(start, end):
   """
   Calculates the end to end distance given two ends of a polymer
   
   Input: start, end coordinates, point array
   -----
   Returns = float(end_to_end_distance)

   """
   
   r = end-start
   return(np.sqrt(np.sum(r**2))) 


def principal_axes(coord):
   """
   Computes the principal axes from moment of inertia.
   
   e1, e2, e3 
   Eigenvectors are sorted by eigenvalues. First eigenvector
   will be the one having higher eigenvalue and so on.

   Input: coord: Input coordinates (nx3) array
   -----

   Returns
   -------
   axis_vectors: 3 x 3 array as [[a1, a2, a3]
                                 [b1, b2, b3]
                                 [c1, c2, c3]]
  
   So, the principal axes are [a1, b1, c1] .. and so on

   However, here the transposed will be returned. 
 
   so, axis = principal_axes(coord)
       1 st axes = axes[0] and so on..
   
   Eigen Vectors are normalized.
   NOTE: Old versions of Gromacs wrote the output data in
         a strange transposed way. As of Gromacs-5.0, the
         output file paxis1.dat contains the x/y/z components
         of the first (major) principal axis for each frame,
         and similarly for the middle and minor axes in 
         paxis2.dat and paxis3.dat. 
 
   """ 
   inertiaTensor = moment_of_inertia(coord)  

   # get eigen values and eigenvector diagonalizing inertia tensor 
   e_val, e_vec = np.linalg.eig(inertiaTensor)

   # get indicies in decreasing order of e_val
   indices = np.argsort(e_val)[::-1]

   # e_vec[:,i] is the eigen vector corresponding to the eigen value e_val[i]
   # Return transposed eigen vector. ie [[--p1--] --> 1st principal axis
   #                                     [--p2--] --> 2nd principal axis
#   return e_vec[:, indices].T         #  [--p3--]]--> 3rd principal axis
   return  e_vec


def anisotropy(coord, **kwargs):
   """
   Measures the shape anisotropy parameter defined as:
 
             k^2 = 1 - 3*((el*e2 + e2*e3 + e3*e1)/(e1+e2+e3)**2)
            
   where e1, e2 and e3 are the eigen values correspondind to 
   the principal axes with e3 >= e2 >= e1.
   It reflects the symmetry and dimensionality of a polymer 
   conformation and lies between values 0 and 1. It reaches 1 for 
   an ideal linear chain and drops to zero for highly symmetric
   conformations. For planar symmetric objects, the relative shape
   anisotropy converges to the value of 1/4.

   asphericity: measures the deviation from the spherical symmetry
               = e3 - 0.5*(e1+e2) ; e3>e2>e1

   For details, follow:  W Janke et al(doi: 10.1063/1.4788616)
   
   Parameters: coord: nd array of shape (nx3). Assume equal mass.
   -------

   Returns: float , value of anisotropy parameter 

   """

   com = np.mean(coord, axis = 0)

   # mass weighted
   mw = 'mw' in kwargs.keys()
   if mw:
      mi = kwargs['mw']
   else: mi = [1]*len(coord)
   total_mass = sum(mi)

  #translate the coordinates to body center of mass
   pos = coord-com

   gyrationTensor = np.zeros((3, 3), dtype=np.float64)

   gyrationTensor[0][0] = (mi*(pos[:, 0] **2 )).sum()
   # xy & yx
   gyrationTensor[0][1] = gyrationTensor[1][0] = (mi*(pos[:, 0] * pos[:, 1])).sum()
   # xz & zx
   gyrationTensor[0][2] = gyrationTensor[2][0] = (mi* (pos[:, 0] * pos[:, 2])).sum()
   # yy
   gyrationTensor[1][1] = (mi*(pos[:, 1]**2 )).sum()
   # yz + zy
   gyrationTensor[1][2] = gyrationTensor[2][1] = (mi*(pos[:, 1] * pos[:, 2])).sum()
   # zz
   gyrationTensor[2][2] = (mi*(pos[:, 2]**2)).sum()

   # divide by total mass
   gyrationTensor /= total_mass
   
   eig_vals = np.linalg.eigvalsh(gyrationTensor)
   # returned eigen values are in ascending order, note that 
   # the eigen values are squared. Gromacs gives sqrt of these values.

   e1, e2, e3 = eig_vals  # e3>e2>e1
#   print(eig_vals) 
   rg_sq = eig_vals.sum()
   rg = np.sqrt(rg_sq)
   anisotropy = 1.0 -3.0 *((e1*e2 + e2*e3 + e3*e1)/(e1+e2+e3)**2)
   asphericity = e3 -0.5*(e1+e2)
   
   return (rg, anisotropy, asphericity)

 
def allign_structure():
   """
   """
   pass 



