#! /usr/bin/env python
"""
   Description: Calculates the radial distribution function.
   *****
   Usage: ./rdf.py -pairs "name O2" "name O2"
   *****
"""
__author = "Hari Sharma <hari.sharma@uconn.edu>"


import math
from math import sqrt, ceil, floor, pow, pi
from os import remove, system, path
import os
from sys import argv, stdout
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from numpy import array, sum
import numpy as np
from numpy import linalg as LA
import structure, topology, mass_data
import timeit

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()

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

def main():
    if(len(argv)) < 2:
      print(parser.print_help()) 
      exit(1)

    binSize = arg.binSize   #default 0.1
    maxrange = arg.maxrange 
    nbins = int(maxrange/binSize)
    maxbin = nbins
    hist   = [0]*(maxbin+1)
    binVolume = [0]*(maxbin+1)
    rdf = []

    nFrame = 0
    nAtomPairs= 0
    totVolume = 0
    process = open('process-rdf.dat', 'w')
    while True:
       grofile = "frame_dump_" + str(nFrame) + ".gro"
       if not path.isfile(grofile) or path.getsize(grofile)==0:
#          print("\n \t File not found. ")
          break
#       elif path.getsize(grofile)==0:
#          print("\n \t File is empty.\n")
#          exit(1)

       molecule = structure.Molecule(filename=grofile)
       selection_string0 = arg.pairs[0]
       selection_string1 = arg.pairs[1]
       if selection_string0==selection_string1:
           coord = molecule.select_atoms(selection_string0)
       box = molecule.box
       nAtoms =  int((len(coord)*(len(coord)-1))/2)  # number of pairs for same types of atoms
       for i in range(len(coord)):
           for j in range(i+1, len(coord)):     
               ri = np.array(coord[i]) 
               rj = np.array(coord[j])
               rij = _vec_pbc(ri, rj, box)
               r = sqrt(np.dot(rij, rij)) 

               bin = int(floor(r/binSize))  # value equal to right edge of bin fall in upper bin.
               if bin <= maxbin:
                  hist[bin] += 1                    
       
       nFrame += 1
       nAtomPairs += nAtoms
       totVolume += box[0]*box[1]*box[2]
       process.write("Processed frame {} \n".format(nFrame))
       process.flush()

    # normalize
    rho = nAtomPairs/totVolume  # average number (pair??) density(N/V)
    for i in range(nbins):
        r1 = i*binSize
        r2 = (i+1)*binSize
        binVolume[i] = (4/3.0)*pi*(pow(r2, 3) -pow(r1,3))     
        r = (i+0.5)*binSize
#        shellVolume = (4./3)*pi*r*r*r
#        binVolume[i] += shellVolume-previousVolume
        hist[i] /= rho*binVolume[i]*nFrame
        rdf.append([r, hist[i]])
 
    rdf = np.array(rdf)
    rdf_file = open('gofr.xvg','w')
    rdf_file.write("@    title \"Radial Distribution Function\" \n")
    rdf_file.write("@    xaxis  label \"r(nm)\" \n")
    rdf_file.write("@    yaxis  label \" freq.\"\n")
    rdf_file.write("@TYPE xy""\n")
    for i in range(len(rdf)):
        rdf_file.write("%8.3f   %8.3f \n" %( rdf[i][0], rdf[i][1])) 

    return rdf    

def args():
    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-bs","--binSize",default = 0.002, help="""Size of bin.Default value is 0.002 nm, just to make it similar
                        to gromacs""",type=float)
    parser.add_argument('-r','--maxrange', default = 35.0, help=""" Maximum range for r. Default value = 35 nm """, type=float)
    parser.add_argument('-f','--grofile', help=""" Starting time for analysis.\n """, type=str)
    parser.add_argument('-pairs', help=""" Need to provide the pairs of atoms for which the rdf has to be calculated
                         .\n """, nargs = 2, type=str)


    return parser, parser.parse_args()

if __name__ =='__main__':
   parser, arg = args()
   main()
