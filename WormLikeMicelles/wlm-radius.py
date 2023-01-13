#! /usr/bin/env python
"""
   This module calculates the radius of an infinite 
   worm like micelle (WLM) slicing the micelle perpendicular
   to the micellar axis. Distance between the center of mass
   of sliced section and head group of micelle gives the radius
   averaged over multiple slice per frame and over trajectory.

   usage: wlm-radius.py -f frame_dump_ -a "name SPh"
"""


from math import sqrt, log
from os import remove, system, path
from sys import argv, stdout
from statistics import stdev
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
from numpy import array, sum, linalg as LA
import os
import timeit

import structure, topology

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()


__author__ = "Hari Sharma"
__email__  = "hari.sharma@uconn.edu"
__version__= "1.0"

class Slice():
    """
       Divide the box into 'n' slices perpendicular to the 
       "direction".
    """
    def __init__(self, coord, box_min, box_max, direction, width=None):
       """
       parameters:
       ----------
       """     

       if width == None:
          self.width = 1.0
       else:
          self.width = width
       self.coord = coord
       self.box_min = box_min
       self.box_max = box_max
       self.box_dim = self.box_max -self.box_min
       self.direction = direction
       self.dir_i   = None
 
       if self.direction.lower() == "x":
               self.dir_i = 0
       elif self.direction.lower() == "y":
               self.dir_i = 1
       else:
               self.dir_i = 2 
       self.nSlices = int(self.box_dim[self.dir_i]/self.width) 
       self.invdelta = self.nSlices/self.box_dim[self.dir_i]

 
    def gen_slice(self):

        slice = [[] for i in range(self.nSlices)]
        for i, pos in enumerate(self.coord):
            r = int((pos[self.dir_i] - self.box_min[self.dir_i])*self.invdelta)%self.nSlices
            slice[r].append(i)
        return slice 

def _dist_sq_pbc(a, b, box):
    """
    Calculates the distance between two points applying pbc.
    a and b are array of single coordinate, size (1x3)
    box is an array with length of box in each dimension
    """
   # applying pbc for a vector between two points
    vec = b -a
    dimension = len(box)
    for i in range(dimension):
        vec[i]= vec[i] - np.rint(vec[i]/box[i])*box[i]
    return np.dot(vec, vec)


def get_com_slices(slices, coord):
    """
    given the index of coordinates, calculates the center of mass.
    parameters:
    ----------
    slices: list of list with indices of atoms within slices
    coord : nx3 array of coordinates  

    returns: array of center of mass coordinates 
    -------
    """
    com = []
    for i in range(len(slices)):
        coord_slice = []
        for j in range(len(slices[i])):
            coord_slice.append(coord[slices[i][j]])
        com.append(np.mean(np.array(coord_slice), axis = 0))
    return np.array(com)

def get_coord_slices(slices, coord):
    """
    given the index of atoms in slices, get the corresponding coordinates
    """
    coord_slice = [[] for i in range(len(slices))]
    for i in range(len(slices)):
        for j in range(len(slices[i])):
            coord_slice[i].append(coord[slices[i][j]])
    return np.array(coord_slice)

def get_radius(com, slices, coord, box):

    if len(com) != len(slices):
       print(" you have unequal number of slices . sth is not \
               quite right")
    distance = []
    for i in range(len(slices)):
        for j in range(len(slices[i])):

            rsq = _dist_sq_pbc(coord[slices[i][j]], com[i], box)
            distance.append(np.sqrt(rsq))
    dist_avg = np.mean(distance)
    return dist_avg 


def get_radius_projection(slices, coord, box, direction):

    """
      calculates the radius using the projection of coordinates in
      each slice to the plane perpendicular to the direction of slicing
      For example, if box is sliced perpendicular to z direction, 
      projection plane is xy. 

      parameters:
      ---------
      slices: index of atoms in each slices, list of lists
      coord : array of coordinates(all coordinates)
      box : array, 1x3:
      direction = int, 0 for x, 1 for y and 2 for z direction
   
    """
    threeD = [0, 1, 2]
    twoD = [x for x in threeD if x != direction]
    distance = []
    box_new = [box[x] for x in twoD]

    for i in range(len(slices)):
        coord_slices = []
        for j in range(len(slices[i])):
            coord_slices.append(coord[slices[i][j]])
        coord_slices = np.array(coord_slices)
        if len(coord_slices) != 0:
           coord_slices_projected = coord_slices[:,[twoD[0], twoD[1]]]
        com = np.mean(coord_slices_projected, axis=0)
        for k in range(len(coord_slices_projected)):
            r = coord_slices_projected[k]- com
            rsq = np.dot(r,r)
#            rsq = _dist_sq_pbc(coord_slices_projected[k], com, box_new)
            distance.append(np.sqrt(rsq))
    dist_avg = np.mean(distance)
    return dist_avg


def main():
    if len(argv) < 2:
       print(parser.print_help())
       exit(1)

    direction = arg.dir
    if direction.lower() == "x":
               dir_i = 0
    elif direction.lower() == "y":
               dir_i = 1
    else:
               dir_i = 2
    width     = arg.width

    
    if arg.gro_file:
       grofile = arg.gro_file 
       
       molecule = structure.Molecule(filename=grofile)
       box_size = molecule.box
       selection_string = arg.atom 
       coordinates_head = molecule.select_atoms(selection_string)
    
       box_min = np.min(coordinates_head, axis = 0)
       box_max = np.max(coordinates_head, axis = 0)

       sl_head = Slice(coordinates_head, box_min, box_max, direction, width)
       slices_head = sl_head.gen_slice()
 
       radius = get_radius_projection(slices_head, coordinates_head, box_size, dir_i)
       print("Radius of WLM: %3.3f"%radius)

   
    elif arg.multi_file:
       radius_all_frames = []

       count = 0
       while os.path.exists(arg.multi_file+"%s.gro"%count):
           grofile = arg.multi_file +str(count) + ".gro"
           print(grofile)
           molecule = structure.Molecule(filename=grofile)
           box_size = molecule.box
           selection_string = arg.atom
           coordinates_head = molecule.select_atoms(selection_string)

           box_min = np.min(coordinates_head, axis = 0)
           box_max = np.max(coordinates_head, axis = 0)

           sl_head = Slice(coordinates_head, box_min, box_max, direction, width)
           slices_head = sl_head.gen_slice()

           radius = get_radius_projection(slices_head, coordinates_head, box_size, dir_i)
           radius_all_frames.append(radius)
           count += 1

    if arg.multi_file: 
     # write the radius of each frame to a file
       rad_wlm = open('radius.xvg','w')
       rad_wlm.write("@    title \"Radius vs time (t)\" \n")
       rad_wlm.write("@    xaxis  label \" t \" \n")
       rad_wlm.write("@    yaxis  label \" radius (nm) \"\n")
       rad_wlm.write("@TYPE xy""\n")      
       for i in range(len(radius_all_frames)):
           rad_wlm.write('{:d} \t {:5.3f} \n'.format(i+1, radius_all_frames[i]))
       stdout.write("="*len("Average radius of WLM = {:5.3f} nm")+"\n")
       stdout.write("Average radius of WLM = {:5.3f} nm \n".format(np.mean(radius_all_frames)))
       stdout.write("="*len("Average radius of WLM = {:5.3f} nm")+"\n")

def args():
    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('-g','--gro_file', help=""" single .gro file for analysis. """, type=str)
    parser.add_argument('-f','--multi_file', help=""" multiple .gro file for analysis without the number and extension. """, type=str)
    parser.add_argument('-a','--atom',     help=""" name of headgroup of wormlike micelles """, type=str)
    parser.add_argument('-d','--dir', dest="dir",    default="z",
                        help=""" direction to slice the box. DEFAULT direction is 'z'. """, type=str)
    parser.add_argument('-w','--width', dest="width",  default = 3.0,  
                        help=""" width of slice in nm. DEFAULT value is '3.0' nm. """, type=float)
    return parser, parser.parse_args()

if __name__ =='__main__':
   parser, arg = args()
   main()


