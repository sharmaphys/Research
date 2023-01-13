#! /usr/bin/env python

"""
#########################################################################
# Written by Hari Sharma 10/24/2017 
# 
# This script calculates order parameter of lipids, and its tilt modulous
# in vesicle using the lipid orientation(tail to head direction) and the 
# normal to each lipids. 
# 
# #######################################################################
"""

from math import sqrt
import math 
from os import remove, system, path
from sys import argv, stdout
import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from numpy import array, sum
import numpy as np
from numpy import linalg as LA
import os
import timeit
from itertools import chain

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()


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

    command1 = "echo %s | trjconv_465_mpi -s %s -f %s -b %s -e %s -n -sep -pbc whole -o frame_dump_.gro > /dev/null"  %(n, tpr, xtc, t1, t2)
    subprocess.call(command1, shell=True)
    return None


def lipid_type(lipidType):
    phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
    phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
    # PCs
    if   lipidType == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
    elif lipidType == "DHPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
    elif lipidType == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
    elif lipidType == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
    elif lipidType == "DEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
    elif lipidType == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
    elif lipidType == "DSPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
    elif lipidType == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
    # PEs
    elif lipidType == "DHPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
    elif lipidType == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
    elif lipidType == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
    elif lipidType == "DSPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
    elif lipidType == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
    elif lipidType == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
    # PPCS
    elif lipidType == "PPCS": bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 AM1-C1A GL2-D1B C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B\n"
    # output legend


    return bond_names


def get_atom_names(bondNames):
    """
    given a bond_names, separates names of atoms in lipid
    and the index of bonded atoms.
    """
    flat = []
    bonds=[x for x in (bondNames.split())]
    atoms=[x.split('-') for x in bonds]
    if bondNames ==[]:
       return flat
    else: 
       for i in atoms:
           for j in i:
              flat.append(j)
    # remove the duplicate names preserving order
    atomList = []
    for i in flat:
        if i not in atomList:
           atomList.append(i)
    # get the index of bonded atoms
    bondedIndex = []
    for i in range(len(bonds)):
        index1st = atomList.index(atoms[i][0])
        index2nd = atomList.index(atoms[i][1])
        bondedIndex.append([index1st, index2nd])
    return(atomList, bondedIndex)
       
def read_gro(grofile, lipidType, surfHead):
   """
   Reads the .gro file and get all the coordinates, coordinates
   of head atoms, number of Lipids, box size and center of vesicle.
   """
   bondNames      = lipid_type(lipidType)
   (atomNames,bondedAtomsIndex) = get_atom_names(bondNames)
   nAtomsPerLipid = len(atomNames)
   nTotalAtoms    = 0
   x1,y1,z1 = 0.0, 0.0, 0.0
   x2,y2,z2 = 0.0, 0.0, 0.0
   atomCoordinates = [[] for i in range(nAtomsPerLipid) ]
   allCoordinates  = []
   headCoordinates = []
   lineCount = 0
   totalLipidAtoms = 0

   for line in open(grofile):
     if lineCount == 1:
      nTotalAtoms = int(line)
     elif lineCount > 1 and lineCount < nTotalAtoms + 2:
        x1 = float(line[20:28])
        y1 = float(line[28:36])
        z1 = float(line[36:44])
        allCoordinates.append([x1,y1,z1])
        if line[5:10].strip() == lipidType and  line[10:15].strip() == surfHead:
           headCoordinates.append([float(line[20:28]),float(line[28:36]), float(line[36:44])])

        for i in range(nAtomsPerLipid):
           if line[5:10].strip() == lipidType and line[10:15].strip()== atomNames[i]:
              x2 = float(line[20:28])
              y2 = float(line[28:36])
              z2 = float(line[36:44])
              atomCoordinates[i].append([x2,y2,z2])
              totalLipidAtoms += 1
     elif lineCount == nTotalAtoms+2:
              box= array(line.split())
     lineCount += 1

   nLipids = int(totalLipidAtoms/nAtomsPerLipid)

   namedCoordinates = [[] for i in range(nAtomsPerLipid)]
   for i in range(nAtomsPerLipid):
       namedCoordinates[i] = array(atomCoordinates[i])
   headCoordinates = array(headCoordinates)
   allCoordinates=array(allCoordinates)
   centerOfMass = [(sum(allCoordinates,axis=0)/float(len(allCoordinates)))]
   com = array(centerOfMass)
   return(namedCoordinates, bondedAtomsIndex, allCoordinates, headCoordinates, com, box, nLipids)

def calc_order_parameter(residueNames, bondIndex, nLipids):
    """
       Calculates the order parameter of each bonded residues for each lipid 
       with respect to the normal.
    """
    nBonds = len(bondIndex)
    normal = [0.0, 0.0, 1.0] 
    P2 = [[] for i in range(nBonds)]
    x, y, z, = 0.0, 0.0, 0.0
    for i in range(nBonds): 
       orderParameter = 0.0
       for j in range(nLipids):
            index1 = bondIndex[i][0]
            index2 = bondIndex[i][1]
            # get the vector between two beads in lipids
            [x,y,z] = residueNames[index1][j]-residueNames[index2][j]
            norm = np.sqrt(x**2 + y**2+z**2)
            normalizedVector = array([x,y,z])/float(norm)
            projection = normal[0]*normalizedVector[0] + normal[1]*normalizedVector[1] + normal[2]*normalizedVector[2]
            orderParameter += projection**2
       P2[i].append(0.5*(3.0*(orderParameter/nLipids) -1.0))            
    return P2

def get_lipid_position(namedCoordinates, headCoordinates, com):
    """ Identifies the position of lipids between inner and outer layer"""
     
    innerLipids, outerLipids = [], []
    dist = []
    midPlaneDist = 0
    nLipids = len(headCoordinates)
    for i in range(nLipids):
        dist.append(np.sqrt(sum(headCoordinates[i] -com)**2))
    midPlaneDist = sum(dist)/float(nLipids)    
   
    for i in range(nLipids):
        if dist[i] <= midPlaneDist:
           innerLipids.append(i)
        else: 
           outerLipids.append(i)        
    return(innerLipids, outerLipids)

def tilt_angle(namedCoordinates,headCoordinates, innerLipids, outerLipids, com, nLipids):
    """ Calculate the tilt angle""" 
    normalInner= []
    normalOuter= []
    tiltAngle  = []

    for i in innerLipids:
         vec = com - headCoordinates[i]
         vecNorm = np.sqrt(sum(vec**2))
         unitVec = vec/float(vecNorm)
         normalInner.append(unitVec)


    for i in outerLipids:
         vec = headCoordinates[i]-com
         vecNorm = np.sqrt(sum(vec**2))
         unitVec = vec/float(vecNorm)
         normalOuter.append(unitVec)

    # get the average of two tail beads of each lipids
    tailCombined = []
    # note that 8 is for C4A, 11 is for C4B, need to automate this\
    tailCombined = (np.add(namedCoordinates[8], namedCoordinates[11]))*0.5
    lipidVecInner = []
    lipidVecOuter = []
    for i in innerLipids:
         vec =  headCoordinates[i]-tailCombined[i]
         vecNorm = np.sqrt(sum(vec**2))
         unitVec = vec/float(vecNorm)
         lipidVecInner.append(unitVec)
    for i in outerLipids:
         vec =  headCoordinates[i]-tailCombined[i]
         vecNorm = np.sqrt(sum(vec**2))
         unitVec = vec/float(vecNorm)
         lipidVecOuter.append(unitVec)


    # get the angle between two unit vectors
    for i in range(len(normalInner)):
        angle = math.acos(np.inner(normalInner[i], lipidVecInner[i]))
        tiltAngle.append(angle)
    for i in range(len(normalOuter)):
        angle = math.degrees(math.acos(np.inner(normalOuter[i], lipidVecOuter[i])))
        tiltAngle.append(angle)

    with open("angle-dist.dat",'w') as f:
        for i in range(len(tiltAngle)):
           f.write("%8.3f \n " %tiltAngle[i])
        
    return(tiltAngle)
     


def main():
    if len(argv) < 2:
       print(parser.print_help())
       exit(1)

    lipidType = "DPPC"
    surfHead  = "NC3"

    if arg.gro_file:
       grofile = arg.gro_file
       (namedCoordinates, bondedAtomsIndex, allCoordinates, headCoordinates, com, box, nLipids) = read_gro(grofile, lipidType, surfHead)
       stdout.write("\n \t Hold on a sec. Working on it.... \n") 
       P2 = calc_order_parameter(namedCoordinates,bondedAtomsIndex,nLipids)  
       (innerLipids, outerLipids)= get_lipid_position(namedCoordinates, headCoordinates, com) 
       tiltAngle = tilt_angle(namedCoordinates,headCoordinates, innerLipids, outerLipids, com, nLipids)

    elif arg.traj_file:
       traj_file() 
       nFrame = 0
       output = open('acf-skip-test.xvg','w')
       output.write("@    title \"Correlation Function\" \n")
       output.write("@    xaxis  label \"Contour Length (s)(nm)\" \n")
       output.write("@    yaxis  label \" C = <u(0).u(s)> \"\n")
       output.write("@TYPE xy""\n") 
       while True:
         grofile = "frame_dump_" + str(nFrame) + ".gro"
         if not path.isfile(grofile) or path.getsize(grofile) == 0:
            break
 
         coordinates = read_gro(grofile)[0]
         if arg.skip_atoms:
            print("\n \t Intermediate atoms are being skipped")
            (coordinates, N) = skip_atoms(coordinates)
            print ("\n \t Number of atoms: %d" %N)
         (CL, CF) =  correlation_function(coordinates)
         print (nFrame)
         remove(grofile)
         nFrame += 1

     #  print(len(contourLength))
       for i in range(len(contourLength)):
          contourLength[i] = sum(contourLength[i])/float(len(contourLength[i]))
          correlationFunction[i] = sum(correlationFunction[i])/float(len(correlationFunction[i]))
          output.write("%8.3f  %8.3f\n" %(contourLength[i], correlationFunction[i])) 


       lp = fit_exponential_decay(contourLength, correlationFunction)
       print ("\t Persistence Length using exp. fit = %8.3f nm \n" %lp)

    end_time = timeit.default_timer()
    totalTimeTaken = end_time-start_time
    stdout.write("\t Time Taken = %i Seconds.\n \n" %totalTimeTaken)



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

    parser.add_argument('-s', '--topol_tpr',
                        help="""
                        topol.tpr file.
                        """,
                        type=str)
                                             
    parser.add_argument('-o','--out_gro',
                        help="""
                        Outputs the coordinates of all micelles.
                        By default outputs only for largest size micelle.
                        """,
                        dest='out_gro', action='store_true')
    
    parser.add_argument('-pbc', '--periodic_boundary',
			help="""
                        Applies the periodic boundary condition for atoms.
			Bt default, it does not apply PBC.
                        """,
                        dest='pbc', action='store_true')

    return parser, parser.parse_args()


if __name__ =='__main__':
   parser, arg = args()
   main()


