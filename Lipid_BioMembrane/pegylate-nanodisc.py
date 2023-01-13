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
    global lipidType,lipidType2, lipidHead, nLipidAtoms,nLipidAtomsType2,  pelResname, pelHead, totalLipids, xPercent
    global resDummy, atomDummy, pelRepeatUnit,pelAtomname, skipnLipid, nPelRequired, resType3

    lipidType = "DPPC"
    lipidHead = "NC3"
    nLipidAtoms = 12
    lipidType2  = "DHPC"
    nLipidAtomsType2 = 8
    resType3    = "GLD"
    nAtomsResType3= 72
    pelResname = "PEL"
    pelHead = "NH3"
    pelAtomname = "O"
    pelRepeatUnit = 90
    resDummy = "DUMM"
    atomDummy = "OOO"
    xPercent = 5
    if len(argv) < 2:
       print(" need to input file name")
       exit(0)
    else:
        grofile = argv[1]
        (nPEL, nLipidsType1, nLipidsType2, nResType3) = gen_peg(grofile)
        topology(nPEL, nLipidsType1, nLipidsType2,nResType3)
    return()



def read_gro(grofile):
   nSurf = 0
   number_of_particles = 0
   lipidHead = "NC3"

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
      if line[5:10].strip() == lipidType and line[10:15].strip() == lipidHead:
         coordinate_head.append([float(line[20:28]),float(line[28:36]), float(line[36:44])])
     line_count += 1

   data = array(coordinate_head)
   array(coordinates)
   center_of_mass = [(sum(data,axis=0)/float(len(coordinate_head)))]
   com = array(center_of_mass)
   coordinate_head = array(coordinate_head)  # head bead of only lipidType
   totalLipids = len(coordinate_head)
   return(coordinates, coordinate_head, com)

def head_group(coordinates, coordinate_head, com):
    """
    Get the index of upper and lower lipids to Pegylate.
    """

    n = len(coordinate_head)
    # generate list of upper and lower lipid index.
    lipidUpper = []
    lipidLower = [] 
    for i in range(n):
       if coordinate_head[i][2]>= com[0][2]+1.0: #trying to avoid for lipids near rim
          lipidUpper.append(i)
       elif coordinate_head[i][2]<= com[0][2]-1.0:
          lipidLower.append(i)   

    nUpper = len(lipidUpper)
    nLower = len(lipidLower)
    print( "\n \t Number of Lipids in upper and lower leaflet that can be pegylated =", nUpper, nLower)

############################################################################
#   In order to get the correct percentage of pegylation, 
#   it is carried out using two steps(halfs)

    skippedHeadAtomUpper = []
    skippedHeadAtomLower = []
    totalLipids = n  # number of total lipids of type lipidType

    nPelRequired = (0.01*xPercent*totalLipids)
    nPelRequiredUpper = int(nPelRequired/2)
    nPelRequiredLower = nPelRequiredUpper
    print(" \t Pegylation required for upper and lower leaflet  = ",nPelRequiredUpper,nPelRequiredUpper)
    skipnLipid_1st_Upper = int(nUpper/float(nPelRequiredUpper))

    nCurrentPelUpper = len(skippedHeadAtomUpper)

    # to pegylate upper leaflet
    for i in range(nUpper):
        if skipnLipid_1st_Upper == 0:
           skippedHeadAtomUpper = lipidUpper  # everything skipped
        elif i%skipnLipid_1st_Upper == 0:
            nCurrentPelUpper += 1
            if nCurrentPelUpper <= nPelRequiredUpper:
               skippedHeadAtomUpper.append(lipidUpper[i])       
    print(" \t Pegylating %d %s Upper Lipids = " %(len(skippedHeadAtomUpper), lipidType))
    
# to pegylate lower leaflet lipids 
   
    skipnLipid_1st_Lower = int(nLower/float(nPelRequiredLower))
    nCurrentPelLower = len(skippedHeadAtomLower)

    for i in range(nLower):
        if skipnLipid_1st_Lower == 0:
           skippedHeadAtomLower = lipidLower  # everything skipped
        elif i%skipnLipid_1st_Lower ==0:
           nCurrentPelLower += 1
           if nCurrentPelLower <= nPelRequiredLower:
              skippedHeadAtomLower.append(lipidLower[i])
    print(" \t Pegylating %s %s Lower Lipids" %(len(skippedHeadAtomLower), lipidType))
    return (skippedHeadAtomUpper, skippedHeadAtomLower)
    
#############################################################################

def gen_peg(grofile):
    with open(grofile) as f1:
        with open("pegylated-nanodisc-"+ str(xPercent) +"%Conc-" + str(pelRepeatUnit) + "repeat.gro", 'w' ) as f2:
          f2.write("PEGylated lipid system\n")
          f2.write("XXXXXXXX \n")
          lines = f1.readlines()
          number_particles = int(lines[1])
          lipid_count = 0
          pel_count = 0 
          atom_count = 0
          (coordinates, coordinate_head, com) = read_gro(grofile)
          (skipped_upper_lipid, skipped_lower_lipid) = head_group(coordinates,coordinate_head,com)

          ithUpperLipid = -1
          ithLowerLipid = -1
          ## pegylate upper lipid
          pos_salt= []
          for i, line in enumerate(lines):
            if line[5:10].strip()==lipidType and line[10:15].strip()==lipidHead:
               ithUpperLipid += 1
               for j in range(len(skipped_upper_lipid)):  
                   if ithUpperLipid == skipped_upper_lipid[j]:
                        ri = array([float(line[20:28]), float(line[28:36]), float(line[36:44])]) # head position
                        f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelHead,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
                        atom_count += 1
                        pel_count += 1
                        pos_x = ri[0]        # for salt coordinates
                        pos_y = ri[1]+0.20   
                        pos_z = ri[2]+0.25
                        pos_salt.append([pos_x, pos_y, pos_z])
                        for k in range(1, nLipidAtoms):
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
                            atom_count += 1
                      
                        r = 0 # z coordinate of head bead
                        for l in range(pelRepeatUnit):
                            r += 0.025
                            ri[2] += r
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,pelAtomname+str(l+1),atom_count+1,ri[0], ri[1], ri[2]))
                            atom_count += 1
             

               ## pegylate  lower lipids
               ithLowerLipid += 1
               for j in range(len(skipped_lower_lipid)):
                   if ithLowerLipid == skipped_lower_lipid[j]:
                        ri = array([float(line[20:28]), float(line[28:36]), float(line[36:44])]) # head position
                        f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelHead,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
                        atom_count += 1
                        pel_count += 1

                        pos_x = ri[0]        # for salt coordinates
                        pos_y = ri[1]+0.20
                        pos_z = ri[2]+0.25
                        pos_salt.append([pos_x, pos_y, pos_z])

                        for k in range(1, nLipidAtoms):
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
                            atom_count += 1

                        r = 0 # z coordinate of head bead
                        for l in range(pelRepeatUnit):
                            r += -0.025
                            ri[2] += r
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,pelAtomname+str(l+1),atom_count+1,ri[0], ri[1], ri[2]))
                            atom_count += 1

	      # for all the rest of lipids not to pegylate, accumulate it
          pegylated = []
          for i in range(len(skipped_upper_lipid)):
              pegylated.append(skipped_upper_lipid[i])
              pegylated.append(skipped_lower_lipid[i])
         
          count = 0
          jj= -1
          resnumber = pel_count
          for i, line in enumerate(lines):
             if i>1 and i <= number_particles+2:
                if line[5:10].strip() == lipidType and line[10:15].strip() == lipidHead: 
                     jj += 1
                     if jj  not in pegylated:
                   
                          x=float(line[20:28])
                          y=float(line[28:36])
                          z=float(line[36:44])
                          f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber+1, lipidType, lipidHead, atom_count+1, x, y, z))
                          atom_count += 1
                          count += 1 
                          resnumber += 1 
                          for k in range(1, nLipidAtoms):
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber,lipidType,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
                            atom_count += 1
                        
          lipid_count  = count
           
          count2 = 0
          ## for different type of lipid
          
          for i, line in enumerate(lines):
             if i>1 and i <= number_particles+2:
                if line[5:10].strip() == lipidType2 and line[10:15].strip() == lipidHead:
                          x=float(line[20:28])
                          y=float(line[28:36])
                          z=float(line[36:44])
                          f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber+1, lipidType2, line[10:15].strip(), atom_count+1, x, y, z))
                          atom_count += 1
                          count2 += 1
                          resnumber += 1
                          for k in range(1, nLipidAtomsType2):
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber,lipidType2,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
                            atom_count += 1

          lipidCountType2 = count2

## for gold, note that all these could be done at once, however to collect all different types of 
## residue it is done separately for each resname type
          count3 = 0
          ## for different type of lipid

          for i, line in enumerate(lines):
             if i>1 and i <= number_particles+2:
                if line[5:10].strip() == resType3:
                          x=float(line[20:28])
                          y=float(line[28:36])
                          z=float(line[36:44])
                          f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber+1, resType3, line[10:15].strip(), atom_count+1, x, y, z))
                          atom_count += 1
                          count3 += 1
                          resnumber += 1
                          for k in range(1, nLipidAtomsType2):
                            f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber,resTYpe3,lines[i+k][10:15], atom_count+1,float(lines[i+k][20:28]), float(lines[i+k][28:36]), float(lines[i+k][36:44])))
                            atom_count += 1

          resCountType3 = count3

          
          pos_salt = array(pos_salt)
          ## generate salt coordinates
          for i in range(len(pos_salt)):
              resnumber += 1
              atom_count += 1
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnumber, "ION", "NA+", atom_count+1, pos_salt[i][0], pos_salt[i][1], pos_salt[i][2]))
              
          f2.write(lines[-1])
    totalLipidPegylated = pel_count
    finalPegPercent = round(totalLipidPegylated/float(totalLipidPegylated + lipid_count), 2)*100
    command = " sed -i 's/XXXXXXXX/%d/g' pegylated-nanodisc-%d%%Conc-%drepeat.gro" %(atom_count,xPercent,pelRepeatUnit)
    subprocess.call(command,shell=True)
    stdout.write(" \n")
    stdout.write("\t PEL Count:- %d \n" %(totalLipidPegylated))
    stdout.write("\t %s Count:-  %d  \n" %(lipidType,lipid_count))
    stdout.write("\t %s Count:-  %d  \n" %(lipidType2,lipidCountType2))
    stdout.write("\n \t FINAL PEGYLATED LIPID PERCENTAGE:-  %d %% \n \n" %finalPegPercent)
    return(pel_count, lipid_count,lipidCountType2, resCountType3)



def topology(pel_count, lipid_count,lipidCountType2, *resCountType3):
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
%s %d
NA+ %d 
''' %(pelRepeatUnit, pel_count, lipidType, lipid_count, lipidType2, lipidCountType2, pel_count ) 	)
     stdout.write("\t Topology written in     ***pel.top*** \n")
     stdout.write("\t Coordinates written in  ***pegylated-nanodisc.gro***  \n\n")


     return ()

if __name__ =='__main__':
    main()
