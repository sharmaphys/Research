#! /usr/bin/env python

"""
Calculates the structure factor of the configuration using
     S(q) = 
"""

import numpy as np
import math
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import structure, topology
from sys import argv, stdout


def static_structure_factor(coordinates, box_size, max_q=None, binMultiplier=None):
    """
    """

# box dimension
    lx = box_size[0]
    ly = box_size[1]
    lz = box_size[2]
# wave vectors    
    dqx = 2*math.pi/float(lx)
    dqy = 2*math.pi/float(ly)
    dqz = 2*math.pi/float(lz)

# max component of q vectors
    if max_q == None:
       max_q = 20.0
    maxQx = math.ceil(max_q/dqx)
    maxQy = math.ceil(max_q/dqy)
    maxQz = math.ceil(max_q/dqz)

# magnitude of maximum q vector,
    maxQ  = math.sqrt(maxQx**2 + maxQy**2 + maxQz**2)
    if binMultiplier==None:
       binMultiplier = 1

    binSize = binMultiplier*min([dqx, dqy, dqz])
    nbins = math.ceil(maxQ/float(binSize))
    
    print(binSize)
    inv_nparticles = 1./float(len(coordinates))

    sq_bin = [[] for i in range(nbins)]
    q_bin = [[] for i in range(nbins)]
  
    # loop over different wave vector
    print("\nCalculating Structure Factor")
    for qx in range(maxQx+1):
        #\r is carriage return for dynamical printing
        print("\t S(q) calculation. Done {:.2f}% \n".format(qx/(maxQx+1)*100), end="\r")  
        for qy in range(maxQy+1):
            for qz in range( maxQz+1):
                q  = np.array([qx*dqx, qy*dqy, qz*dqz])
                normQ = np.linalg.norm(q)
                if normQ <= max_q:
                   # determine which bin q lies on
                   ibin = math.floor(normQ/float(binSize))   
                   q_bin[ibin].append(normQ)


                   cos_temp = 0.0
                   sin_temp = 0.0

                   # loop over particles
                   for i in range(len(coordinates)):
                       r = np.array(coordinates[i])
                       qdotR = np.dot(q, r)
                       cos_temp += math.cos(qdotR)
                       sin_temp += math.sin(qdotR)
                   sq_bin[ibin].append(cos_temp*cos_temp + sin_temp*sin_temp)            
    q_bin_nonempty = [x for x in q_bin if x] 
    sq_bin_nonempty = [x for x in sq_bin if x]  

    # get the average values normalized by number of particles
    q_bin_avg  = [sum(x)/float(len(x)) for x in q_bin_nonempty]
#    q_bin_avg  = [i*binSize for i in range(len(sq_bin_nonempty))]
    sq_bin_avg = [(sum(x)/float(len(x)))*inv_nparticles for x in sq_bin_nonempty]

    if len(q_bin_avg) != len(sq_bin_avg):
       raise ValueError('Something is not quite right in S(q) calculation.\
                         Unequal bin length for q and S(q)')
       exit
    return q_bin_avg, sq_bin_avg
   

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
     
