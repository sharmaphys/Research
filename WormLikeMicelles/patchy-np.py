#! /usr/bin/env python3.6
"""
Created on Fri May 18 20:36:32 2018

@author: has12010
"""

import numpy as np
import datetime
import math

def gen_itp():
    pass
def write_gro():
    pass

def fcc(latticeConstant):
    #
    a = latticeConstant
    #basis vectors
    r =  [[0.0,   0.0,   0.0],
         [a*0.5, a*0.5, 0.0],
         [a*0.5, 0.0, a*0.5],
         [0.0, a*0.5, a*0.5]]
    return r

def fibonacci_sphere(n, r, center=None):
    """
       Gives the n points uniformly distributed in a
       sphere of radius r.
          
       inputs: n = integer,  number of points to be generated
               r = float, radius of sphere
               center = array, center of sphere

       Golden angle = pi*(3-sqrt(5)), the most irrational angle
              theta = golden_angle* i
              z[i] = (1 -1./n)*(1 - 2*i/(n -1))
              radius[i] = sqrt(1-(z[i])**2)
    """
    if center == None:
       center = np.array([0.0, 0.0, 0.0 ])
    else:
       center = np.array(center)
    golden_angle = math.pi*(3 - math.sqrt(5))

    theta = np.arange(n)*golden_angle
    
    z = np.linspace(1-1.0/n, 1.0/n -1, n)  # min of z value when i = 0 => 1-1.0/n, 
 				           # max is when i = n-1, => 1.0/n -1
    radius = np.sqrt(1-z*z)

    points = np.zeros((n, 3))
    points[:,0] = radius * np.cos(theta)*r + center[0] # scaled with desired radius r and center
    points[:,1] = radius * np.sin(theta)*r + center[1]
    points[:,2] = z*r + center[2]
  
    return points


def gen_np(radius, beadsPerLigand, latticeConstant=None):
    if latticeConstant ==None:
       a = 0.408  # for gold fcc lattice, 4 atoms per unit cell
    else:
       a = latticeConstant
    atomPos = []
    
    ## basis vectors
    r = [[0.0,   0.0,   0.0],
         [a*0.5, a*0.5, 0.0],
         [a*0.5, 0.0, a*0.5],
         [0.0, a*0.5, a*0.5]]
    
    for i in range(-8,9):
        for j in range(-8,9):
            for k in range(-8,9):
                
                #for each basis vectors
                for l in range(4):
                    xx = i*a +  r[l][0] 
                    # still do not know why there should be (j+0.5)
                    yy = (j+0.5)*a +  r[l][1] 
                    zz = k*a +  r[l][2] 
                    atomPos.append([xx, yy, zz])
    atomPos = np.array(atomPos)                
    npCore = []                
    for i in range(len(atomPos)):
#   to get the truncated octahedron NP, use the first two commented line
#        dist = abs(atomPos[i][0]) + abs(atomPos[i][1]) + abs(atomPos[i][2])
#        if dist <= 1.4143*radius and abs(atomPos[i][0])<=radius and abs(atomPos[i][1])<=radius and abs(atomPos[i][2]) <= radius:
        rsq = atomPos[i][0]**2 + atomPos[i][1]**2 + atomPos[i][2]**2
        if rsq <= radius**2 :# and abs(atomPos[i][0])<=radius and abs(atomPos[i][1])<=radius and abs(atomPos[i][2]) <= radius :
           npCore.append([atomPos[i][0], atomPos[i][1], atomPos[i][2]])
    return np.array(npCore)                


def get_bonds(npCore, latticeConstant):
    n = len(npCore)
    numBonds = np.zeros((n, 1))
    for i in range(n):
        for j in range(i+1, n):
            dist = np.linalg.norm(npCore[j]- npCore[i])   
            # distance between nearest neigbor in fcc = (1/2)*sqrt(2)*lattice constant 
            if dist*dist <= 0.5*(latticeConstant**2):
               numBonds[i] += 1
               numBonds[j] += 1
    shell_atoms = []
    for i in range(n):
        if numBonds[i] != 12:
           shell_atoms.append(npCore[i])
    return shell_atoms 


def get_shell_atoms(npCore):
    """
    Given a nanoparticle core, gets the shell atoms. ie.the 
    atoms on the surface of the nanoparticle.
    """
    pass


def write_gro(atomPos):
    resname = "GLD"
    atomname = "N0"
    box_size = 15.0
    now = datetime.datetime.now().strftime("%m/%d/%Y %H:%M")
    atomPos = np.array(atomPos)
    f=open('sphere.gro','w')
    f.write("Nanoparticle created by Hari on %s \n" %now)
    f.write("%d \n" %len(atomPos))
    for i in range(len(atomPos)):
        f.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(i+1, resname, atomname, i+1, atomPos[i][0], atomPos[i][1], atomPos[i][2]))    
    f.write('%8.3f%8.3f%8.3f\n' %(box_size,box_size,box_size))
    f.close
    
    
def main():
#    atomPos = gen_np(1.,10)
    atomPos = fibonacci_sphere(800, 5) 
    print(len(atomPos))
#    shellAtoms = get_bonds(atomPos, 0.408)
    write_gro(atomPos)
#    write_gro(atomPos)
     
if __name__== "__main__":
   main()                                                              

