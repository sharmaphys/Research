#! /usr/bin/env python

"""
Calculates the structure factor of the configuration using
     S(q) = 
"""

import numpy as np
import math
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import structure
from sys import argv, stdout


def static_structure_factor_2d(coordinates, max_q):
    """
    """

    ete_vector = np.abs(np.array(coordinates[-1]) -np.array(coordinates[0]))
    # sorts in ascending order, gives the indicies of element after sorting.
    print(ete_vector)
    indicies = np.argsort(ete_vector)
    ind1 = indicies[0]
    ind2 = indicies[1]
    ind3 = indicies[2] # index of max element
    print("ind1:", ind1, "ete_ind1", ete_vector[ind1], "ind2:", ind2, "ete_ind2", ete_vector[ind2])
    lz = ete_vector[ind3]
    print("lz",lz)
    # wave vectors    
    dq = 2*math.pi/float(lz)

# max component of q vectors
    if max_q == None:
       max_q = 40.0

# magnitude of maximum q vector,
#    if binMultiplier==None:
#       binMultiplier = 1
    binSize = dq
    nbins = math.ceil(max_q/float(binSize))
#    nbins = 2*n
    
    print(binSize)
    inv_nparticles = 1./float(len(coordinates))

    sx_bin = [[] for i in range(2*nbins+1)]
    sy_bin = [[] for i in range(2*nbins+1)]
#    sz_bin = [[] for i in range(2*nbins+1)]

    q_bin = [[] for i in range(2*nbins+1)]
  
    # loop over different wave vector
    print("\nCalculating Structure Factor")
    for i in range(1, nbins):
#        print("\t S(q) calculation. Done {:.2f}%".format(i/(max_q+1)*100))
        q_i  = i*binSize
        if q_i <= max_q:
          # determine which bin q lies on
          ibin = math.floor(q_i/float(binSize))   
          q_bin[ibin].append(q_i)

          cos_tempX = 0.0
          sin_tempX = 0.0
          cos_tempY = 0.0
          sin_tempY = 0.0

          # loop over particles
          for j in range(len(coordinates)):
             x = np.array(coordinates[j][ind1])
             y = np.array(coordinates[j][ind2])
#             z = np.array(coordinates[j][ind3])

             qX = x*q_i
             qY = y*q_i
#             qZ  = z*q_i

             cos_tempX += math.cos(qX)
             sin_tempX += math.sin(qX)
             cos_tempY += math.cos(qY)
             sin_tempY += math.sin(qY)
#             cos_tempZ += math.cos(qZ)
#             sin_tempZ += math.sin(qZ)

          sx_bin[ibin].append(cos_tempX*cos_tempX + sin_tempX*sin_tempX)
          sy_bin[ibin].append(cos_tempY*cos_tempY + sin_tempY*sin_tempY)                   
 #         sz_bin[ibin].append(cos_tempZ*cos_tempZ + sin_tempZ*sin_tempZ)

    q_bin_nonempty = [x for x in q_bin if x] 
    sx_bin_nonempty = [x for x in sx_bin if x]  
    sy_bin_nonempty = [x for x in sy_bin if x]
#    sz_bin_nonempty = [x for x in sz_bin if x]

    # get the average values normalized by number of particles
    q_bin_avg  = [sum(x)/float(len(x)) for x in q_bin_nonempty]
#    q_bin_avg  = [i*binSize for i in range(len(sq_bin_nonempty))]
    sx_bin_avg = [(sum(x)/float(len(x)))*inv_nparticles for x in sx_bin_nonempty]
    sy_bin_avg = [(sum(x)/float(len(x)))*inv_nparticles for x in sy_bin_nonempty]
#    sz_bin_avg = [(sum(x)/float(len(x)))*inv_nparticles for x in sz_bin_nonempty]


    if len(q_bin_avg) != len(sx_bin_avg):
       raise ValueError('Something is not quite right in S(q) calculation.\
                         Unequal bin length for q and S(q)')
       exit
    return q_bin_avg, sx_bin_avg, sy_bin_avg
   

def main():
    if len(argv) < 2:
       print(parser.print_help())
       exit(1)

    if arg.gro_file:
       nFrame = 1
       grofile = arg.gro_file

       sf = open('sf.xvg','w')
       sf.write("@    title \"Structure Factor\" \n")
       sf.write("@    xaxis  label \" q (1/nm)\" \n")
       sf.write("@    yaxis  label \" S(q) \"\n")
       sf.write("@    yaxis  label \" q*S(q) \"\n")
       sf.write("@TYPE xy""\n")
 
       molecule = structure.Molecule(filename=grofile)
       coordinates = molecule.get_coord_all()
   #    coordinates = molecule.select_atoms("resname DHPC and name NC3")
       box_size = molecule.box
       print(box_size)
       q_bin, sq_bin = static_structure_factor(coordinates, box_size, max_q=None, binMultiplier=None)
 
       # write values excluding q = 0.
       for i in range( len(q_bin)):
           sf.write(("{:.3f}" + "\t" +"{:.3f}" + "\t" + "{:.3f}" + "\n").format(q_bin[i], sq_bin[i], q_bin[i]*sq_bin[i]))


def args():
    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-f","--traj_file",
                        help=""".xtc(trajectory) file for analysis.""",
                        type=str)

    parser.add_argument('-g','--gro_file',
                        help="""
                        Single .gro file for analysis.
                        """,
                        type=str)

    parser.add_argument('-b','--begin_time',
                        help="""
                        Starting time for analysis.
                        \n        
                        """,
                        type=int)

    parser.add_argument('-e','--end_time',
                        help="""
                        End time for analysis.
                        """,
                        type=int)

    parser.add_argument('-n','--index_no',
                        help="""
                        Index number for polymer chain.
                        """,
                        type=int)
    return parser, parser.parse_args()

if __name__ =='__main__':
   parser, arg = args()
   main()
     
