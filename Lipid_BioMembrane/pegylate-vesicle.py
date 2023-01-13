#! /usr/bin/env python
"""

###########################################################################
# Description:  This script pegylates ie replace lipids by PEG lipid 
#               given a vesicle with percentage of required pegylation.
##########################################################################

"""
from math import sqrt,ceil,floor
from os import remove, system, path
from sys import argv, stdout
import subprocess
import numpy as np
from numpy import array, sum
#from termcolor import colored,cprint



def main():
    global lipidType, lipidHead, nLipidAtoms, pelResname, pelHead, totalLipids, xPercent
    global resDummy, atomDummy, pelRepeatUnit,pelAtomname, skipnLipid, nPelRequired

    lipidType = "DPPE"
    lipidHead = "NH3"
    nLipidAtoms = 12
    pelResname = "PEL"
    pelHead = "NH3"
    pelAtomname = "O1"
    pelRepeatUnit = 45
    #skipnLipid = 20  # number of lipids to skip while pegylating.
    resDummy = "DUMM"
    atomDummy = "OOO"
    xPercent = 24
    if len(argv) < 2:
       print " need to input file name"
       exit(0)
    else:
        grofile = argv[1]
        (nPEL, nLipids) = gen_peg(grofile)
        topology(nPEL, nLipids)
    return()



def read_gro(grofile):
   nSurf = 0
   number_of_particles = 0
   lipidHead = "NH3"
   x,y,z = 0.0, 0.0, 0.0
   coordinates = []
   coordinate_head = []
   line_count = 0
   for line in open(grofile):
     if line_count == 1:
      number_of_particles = int(line)
      nSurf = number_of_particles/nLipidAtoms
     elif line_count > 1 and line_count < number_of_particles + 2:
      x = float(line[20:28])
      y = float(line[28:36])
      z = float(line[36:44])
      coordinates.append([x,y,z])
      if line[10:15].strip() == lipidHead:
         coordinate_head.append([float(line[20:28]),float(line[28:36]), float(line[36:44])])
     line_count += 1

   data = array(coordinate_head)
   array(coordinates)
   center_of_mass = [(sum(data,axis=0)/float(len(coordinate_head)))]
   com = array(center_of_mass)
   coordinate_head = array(coordinate_head)
   totalLipids = len(coordinate_head)
   return(coordinates, coordinate_head, com)

def head_group(coordinates, coordinate_head, com):
    """
    Get the distance between the com and head atom of surfactants
    and output the outer spheres of headgroup only.
    """

    nParticles = len(coordinate_head)
    with open('all-head.gro', 'w') as out:
         out.write(" Head groups  \n")
         out.write(" %s \n"%(nParticles+1))  # plus 1 is form dummy atom
         lipid_count = 0
         atom_count = 0
         for i in range(nParticles):
              j = i*nLipidAtoms
              x = coordinates[j][0]
              y = coordinates[j][1]
              z = coordinates[j][2]
              out.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(lipid_count+1,lipidType,lipidHead,atom_count+1, x, y, z))
              lipid_count +=1
              atom_count += 1
         out.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(lipid_count,resDummy,atomDummy,atom_count+1, com[0][0], com[0][1], com[0][2]))
         out.write("  35.00  35.00  35.00 \n")   


    # output only outer sphere of head atoms.
    distance = [] 
    for i in range(nParticles):
       dist = np.sqrt(np.sum((coordinate_head[i]-com)**2))
       distance.append(dist)

    meanDistance = ceil(np.mean(array(distance)))
    print "\n \t Mean distance between head atoms in vesicle =", meanDistance
    headAtomOuter = []
    for i in range(nParticles):
       dist = np.sqrt(np.sum((coordinate_head[i]-com)**2))
       distance.append(dist)
       
       if dist > meanDistance:
           headAtomOuter.append(i)
    headAtomOuter=array(headAtomOuter)
    print "\t Number of Head Atoms in Outer Sphere =", len(headAtomOuter)

    # generate .gro file with outer head beads only
    n = len(headAtomOuter)
    with open('outer-head.gro', 'w') as out:
         out.write(" Outer Head groups  \n")
         out.write(" %s \n"%(n+1))  # plus 1 is form dummy atom
         lipid_count = 0
         atom_count = 0
         for i in range(n):
                j = headAtomOuter[i]*nLipidAtoms
                x = coordinates[j][0]
                y = coordinates[j][1]
                z = coordinates[j][2]
                out.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(lipid_count+1,lipidType,lipidHead,atom_count+1, x, y, z))
                lipid_count +=1
                atom_count += 1
         out.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(lipid_count+1,resDummy,atomDummy,atom_count+1, com[0][0], com[0][1], com[0][2]))
         out.write("   35.000 35.000 35.000 \n")


############################################################################
#   In order to get the correct percentage of pegylation, 
#   it is carried out using two steps(halfs)

    skippedHeadAtomOuter = []
    totalLipids = nParticles  # number of total lipids

    nPelRequired = (0.01*xPercent*totalLipids)
    skipnLipid_1st = floor(n/float(nPelRequired))
    skipnLipid_2nd = ceil(n/float(nPelRequired))
    firstHalf     = int(ceil(n/2))
    for i in range(firstHalf):
        if skipnLipid_1st == 0:
           skippedHeadAtom = headAtomOuter
        elif i%skipnLipid_1st ==0:
           skippedHeadAtomOuter.append([headAtomOuter[i]])       
    
    reqdMore = int(nPelRequired) - len(skippedHeadAtomOuter)

    nCurrentPel = len(skippedHeadAtomOuter)
    print("\n \t Required PEL: %d, Current PEL: %d, Needed : %d more  " %(int(nPelRequired),nCurrentPel,reqdMore))

    for i in range(firstHalf, n):
      #   while(nCurrentPel <= nPelRequired):
            if skipnLipid_2nd ==0:
               skippedHeadAtom = headAtomOuter
            elif i%skipnLipid_2nd ==0:
#              skippedHeadAtomOuter.append([headAtomOuter[i]])
              nCurrentPel += 1
              if nCurrentPel <= nPelRequired:
                  skippedHeadAtomOuter.append([headAtomOuter[i]])
    print(" \t Pegylating %s %s Lipids" %(len(skippedHeadAtomOuter), lipidType))

#############################################################################
    return (array(skippedHeadAtomOuter))

def gen_peg(grofile):
    with open(grofile) as f1:
         with open('peg-vesicle.gro', 'w') as f2:
	      f2.write("PEGylated lipid system\n")
	      f2.write("XXXXXXXX \n")
	      lines = f1.readlines()
	      number_particles = int(lines[1])
	      lipid_count = 0
	      pel_count = 0
	      atom_count = 0
	      (coordinates, coordinate_head, com) = read_gro(grofile)
    #          totalLipids = len(coordinate_head)
	      skipped_outer_sphere = head_group(coordinates,coordinate_head,com)
	      for i, line in enumerate(lines):
	          for j in range(len(skipped_outer_sphere)):
	             if i>1: 
	              if i==skipped_outer_sphere[j]*nLipidAtoms + 2:   # plus 2 as real atom starts from line 2
	                    #print i, skipped_outer_sphere[j]
	                    ri = array([float(line[20:28]), float(line[28:36]), float(line[36:44])])
#	                    rj = array([float(lines[i+ nLipidAtoms-1][20:28]), float(lines[i + nLipidAtoms-1][28:36]), float(lines[i+nLipidAtoms-1][36:44])])
                        rVec =  ri-[com[0][0],com[0][1],com[0][2]]   # vector from center to head

#	                    rVec =  ri-rj   # vector from tail to head
	                    rNorm = np.sqrt(sum(rVec**2))
	                    unitVec = rVec/float(rNorm)
 
#          
	                    f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelHead,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
	                    atom_count += 1
	                    pel_count += 1

	                    for k in range(1, nLipidAtoms):

	                      f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
	                      atom_count += 1
                      
			    r = 0
	                    for l in range(pelRepeatUnit):
	                         r += unitVec*0.03
	                         ri += r
	                         f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,pelAtomname,atom_count+1,ri[0], ri[1], ri[2]))
	                         atom_count += 1


	      # for all the rest of lipids not to pegylate, accumulate it
	      count = 0
	      resnumber = pel_count
	      for i, line in enumerate(lines):
	           if i>1 and i <= number_particles+2:
	               if line[5:10].strip() == lipidType and line[10:15].strip() == lipidHead:
	                  j = (i -2)/float(nLipidAtoms)
	                  if j not in skipped_outer_sphere:
	                      x=float(line[20:28])
	                      y=float(line[28:36])
	                      z=float(line[36:44])
	                      f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber+1, lipidType, pelHead, atom_count+1, x, y, z))
	                      atom_count += 1
	                      count += 1 
	                      resnumber += 1 
	                 #     count1 += 1
	                      for k in range(1, nLipidAtoms):
	                        f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber,lipidType,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
	                        atom_count += 1
                        
	      lipid_count = count  
	      f2.write("35.000 35.000 35.000 \n")
    finalPegPercent = round(len(skipped_outer_sphere)/float(len(skipped_outer_sphere) + lipid_count), 2)*100
    command = " sed -i 's/XXXXXXXX/%d/g' peg-vesicle.gro" %(atom_count)
    subprocess.call(command,shell=True)
    stdout.write(" \n")
    stdout.write("\t PEL  Count:-  %d  \n"%(len(skipped_outer_sphere)))
    stdout.write("\t %s Count:-  %d  \n" %(lipidType,lipid_count))
    stdout.write("\n \t FINAL PEGYLATED LIPID PERCENTAGE:-  %d %% \n \n" %finalPegPercent)
    return(pel_count, lipid_count)



def topology(pel_count, lipid_count):
     with open ('pel.top','w') as top:
          top.write( '''
; Last updated 20 july 2011
;
; Note: when using this topology, please make sure to change the
;       following line in you martini_v2.x.itp file:
;
; SN0   SN0     1       0.85338E-01     0.53946E-03 ; 75almost attractive, s=0.43
; (instead of SN0 SN0 1 0.66375E-01     0.41957E-03 ; 75intermediate, s=0.43)
;
; AND	
;
;  P4   SN0     1       0.17246E-00     0.18590E-02 ; semi attractive
; (instead of  P4 SN0 1 0.15091E-00     0.16267E-02 ; intermediate) 
;
; This is required to enhance the self-interaction between the PEO monomers
; and the interaction between PEO monomers and water
; to the level of an Nda particle, see the reference above and also
; Lee & Larson, J. Phys. Chem. B, 2009, 113 (40), pp 13202-13207
; Interaction with other particles types, initially a SNa type particle was chosen,
; but later a SN0 type particle was shown to behave better, see
; Lee & Pastor, J. Phys. Chem. B, 2011, 115, 7830-7837
;
; Be cautious if you have other particles in your system of type 'SN0',
; their self- and water-interaction will also be affected.

#include "martini_v2.0_for_pegylated.itp"
#include "PEL%d.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"


[ system ]
PEGylated bilayer

[ molecules ]
PEL %d
%s %d
NA+ %d 
''' %(pelRepeatUnit, pel_count, lipidType, lipid_count, pel_count ) 	)
     stdout.write("\t Topology written in     ***pel.top*** \n")
     stdout.write("\t Coordinates written in  ***peg-vesicle.gro***  \n\n")


     return ()

if __name__ =='__main__':
    main()
