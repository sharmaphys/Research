#! /usr/bin/env python
"""
 This script calculates the hydrogen bond between the donar and acceptor atoms
 if it satisfies the geometrical (distance and angle) criteria, ie if the distance
 between the donor and acceptor atom is less than or equal to 0.35 nm and the 
 angle between the hydrogen-donor-acceptor is less than 30 degrees
 Cutoff implies the maximum distance between donor and acceptor atom
 and the maximum deviation from straight 180 degrees angle.

 Overall calculation is based on grid search. Gridding is done based on 
 the dimensions of all coordinates. Positions of acceptors in each grid
 is found and look for donors in all neighbor cells that satisfy both 
 distance and angle cutoff criteria.
"""
__author__ = "Hari Sharma <hari.sharma@uconn.edu Date: 5/3/2018>"

import math
import itertools
from sys import argv,stdout
import numpy as np
import timeit
import structure
#numpy.set_printoptions(threshold=numpy.nan)
start_time = timeit.default_timer()


class Grid():
     def __init__(self, coord, box_min, box_max, cutoff=None):

         if cutoff == None:
            self.cutoff = 0.35
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


def get_hbond(coord_acceptor, coord_donor, coord_donor_all, coord_all):
    """
    Parameters: 
    ----------
         coord_acceptor: coordinates of acceptor
         coord_donor: coordinates of donor
         coord_donor_all: coordinates of all atom of solvent(donor)
    """
    box_min = np.min(coord_all, axis=0)
    box_max =  np.max(coord_all, axis=0)
    box_dim   = box_max-box_min

    gd = Grid(coord_donor, box_min, box_max)
    cellsize = gd.cutoff
    invdelta = gd.invdelta
    grid = gd.gen_grid()
    nGrids = gd.ngrids

    distcutoff = 0.35
    distcutoff_sq  = distcutoff**2
    angleCutoff = 30

    # loop over acceptor
    hbond = []
    # accumulate number of donors(eg. OW of SOL) within cutoff of acceptor
    first_hydration_shell = []
    for i in range(len(coord_acceptor)):

        neigh_atom_ndx = []
        # index of acceptor atom in cell
        cell = np.array([ int((coord_acceptor[i][0]-box_min[0])*invdelta[0])%nGrids[0],  \
                          int((coord_acceptor[i][1]-box_min[1])*invdelta[1])%nGrids[1],  \
                          int((coord_acceptor[i][2]-box_min[2])*invdelta[2])%nGrids[2] ])
        neigh_atom_ndx = get_neighbor(cell, grid, nGrids)
        for j in range(len(neigh_atom_ndx)):
             distsq = _dist_sq_pbc(coord_acceptor[i], coord_donor[neigh_atom_ndx[j]], box_dim)
             # check distance criteria
             if distsq <= distcutoff_sq:
              # if distance criteria is satisfied, check angle criteria
              # angle between hydrogen-donor-acceptor should be <= 30 
              # i. e max deviation from linearity [which is 0] is 30. 


              # get hydrogen atom coordinates bonded to jth heavy atom (OW) of SOL
              # will have indicies j + 1 and j + 2
              first_hydration_shell.append(coord_donor[neigh_atom_ndx[j]])
              for k in range(1,3):

                  hd = _vec_pbc(coord_donor[neigh_atom_ndx[j]], coord_donor_all[neigh_atom_ndx[j]*3 + k], box_dim)
                  ad = _vec_pbc(coord_donor[neigh_atom_ndx[j]],  coord_acceptor[i], box_dim)
                  angle  = _angle_between(hd, ad)
                  if cell[0] == nGrids[0]:
                     print(cell, coord_acceptor[i], distsq, hd, ad, angle)

                  if angle <= angleCutoff:
                     hbond.append([i, j])

    print("\n\t Number of H-bonds: ",len(hbond))
    print("\t Number of donors in 1st Hydration Shell: ",len(first_hydration_shell))
    print("\t Using : ", nGrids[0],'x',nGrids[1],'x',nGrids[2], " number of grids with cell size = %3.2f"%cellsize )
    return hbond

def main():

    molecule = structure.Molecule(filename= argv[1])

    # get coordinate of hbond acceptor (polymer)
    coord_acceptor = molecule.get_coord_atom("O1", resname="OLEA")
    box_min_acceptor =  np.min(coord_acceptor, axis=0)
    box_max_acceptor =  np.max(coord_acceptor, axis=0)

    # get coordinate of donor (from Solvent)
    coord_donor = molecule.get_coord_atom("OW", resname="SOL")
    box_min_donor =  np.min(coord_donor, axis=0)
    box_max_donor =  np.max(coord_donor, axis=0)

    # get coordinate of whole solvent with hydrogen
    coord_sol = molecule.get_coord_resname("SOL")
    coord_all = molecule.get_coord_all()  
#    print(np.min(coord_all, axis=0))
    hbond = get_hbond(coord_acceptor, coord_donor, coord_sol, coord_all)

if __name__=="__main__":
   main()
   end_time = timeit.default_timer()
   print('\t Time taken : %3.3f Seconds \n'%(end_time - start_time))

