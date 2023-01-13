#! /usr/bin/python

###############################################################
#  this script replaces the DPPC lipid by PEGylated lipid (PEL) 
#  and adds peg monomers on the head group of replaced DPPC
###############################################################


from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
from numpy import array, sum
#from termcolor import colored,cprint

skip_dppc = 12  ## to get the correct number of grafting density of peg
beads_dppc = 12
beads_dhpc = 8
lipid_type = "DPPC"
atom1 = "NC3"
atom2 = "C4B"
pelLipidHead = "NH3"
pelResname = "PEL"
pelAtomName = "O1"
repeatUnit = 90       # repeat unit of peg monomer
with open('nanodisc.gro') as f1:
    with open('peg.gro', 'w') as f2:
      f2.write("PEGylated lipid system\n")
      f2.write("XXXXXXXX \n")
      lines = f1.readlines()
      dppc_count = 0
      pel_count = 0
      atom_count = 0
      for i, line in enumerate(lines):
         number_particles = lines[1]
         if line[5:10].strip() == lipid_type and line[10:15].strip() == atom1:
#            NC3_Z = float(line[36:44])
#            C4B_Z = lines[i+11][36:44]
#            print NC3_Z
#            print C4B_Z
            if dppc_count%skip_dppc == 0:
              NC3_Z = float(line[36:44])
              C4B_Z = float(lines[i+11][36:44])
              Z1 = NC3_Z
#              for k in range(repeatUnit):
#                if NC3_Z > C4B_Z:
#                 Z1 += 0.20
#                 f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelAtomName,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), Z1))
#                 atom_count += 1
#            
#                elif NC3_Z < C4B_Z:
#                 Z1 -= 0.20
#                 f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelAtomName,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(Z1)))
#                 atom_count += 1
#           
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+1,pelResname,pelLipidHead,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
              atom_count += 1
              pel_count += 1

              for j in range(1, beads_dppc):

                f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,lines[i+j][10:15], atom_count+1,float(lines[i+j][20:28]), float(lines[i+j][28:36]), float(lines[i+j][36:44])))
                atom_count += 1

              for k in range(repeatUnit):
                if NC3_Z > C4B_Z:
                 Z1 += 0.20
                 f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,pelAtomName,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), Z1))
                 atom_count += 1

                elif NC3_Z < C4B_Z:
                 Z1 -= 0.20
                 f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count,pelResname,pelAtomName,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(Z1)))
                 atom_count += 1

            dppc_count += 1
      count = 0
      count2 = 0
      lipid_count_DPPC = 0
      lipid_count_DHPC = 0
      lipid_type1 = "DPPC"
      lipid_type2 = "DHPC"
      for i, line in enumerate(lines):
         if line[5:10].strip() == lipid_type1 and line[10:15].strip() == atom1:
         #  print lipid_count_DPPC
           if lipid_count_DPPC%skip_dppc != 0:
              count += 1 
#              print lipid_count_DPPC%skip_dppc
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(count+pel_count,lipid_type1,atom1,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
#              lipid_count_DPPC += 1
              atom_count += 1
              for j in range(1, beads_dppc):

                f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(count+pel_count,lipid_type1,lines[i+j][10:15],atom_count+1,float(lines[i+j][20:28]), float(lines[i+j][28:36]), float(lines[i+j][36:44])))
                atom_count += 1
           lipid_count_DPPC += 1        
              #print lipid_count_DPPC 


         elif line[5:10].strip() == lipid_type2 and line[10:15].strip() == atom1:
              count += 1
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+count,lipid_type2,atom1,atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
             # lipid_count_DHPC += 1 
              atom_count += 1

              for j in range(1, beads_dhpc):

                f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(pel_count+count,lipid_type2,lines[i+j][10:15],atom_count+1,float(lines[i+j][20:28]), float(lines[i+j][28:36]), float(lines[i+j][36:44])))
                atom_count += 1
              lipid_count_DHPC += 1
      f2.write("35.000 35.000 35.000 \n")
command = " sed -i 's/XXXXXXXX/%d/g' peg.gro" %(atom_count)
subprocess.call(command,shell=True)
stdout.write(" \n")
stdout.write("  PEL  Count:-  %d  \n"%(pel_count))
stdout.write("  DPPC Count:-  %d  \n"%(count-lipid_count_DHPC))
stdout.write("  DHPC Count:-  %d  \n\n"%(lipid_count_DHPC))
with open ('pel.top','w') as lipidType:
  lipidType.write( '''
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


  #include "martini_v2.0-PEG.itp"
  #include "PEL%d.itp"
  #include "martini_v2.0_lipids.itp"
  #include "martini_v2.0_ions.itp"


  [ system ]
  PEGylated bilayer

  [ molecules ]
  PEL %d
  DPPC %d
  DHPC %d
  DPPC %d
  DHPC %d
  NA+ %d 
  ''' %(repeatUnit, pel_count, (count-lipid_count_DHPC)/2, (lipid_count_DHPC)/2, (count-lipid_count_DHPC)/2,(lipid_count_DHPC)/2,pel_count )

  )
 # lipidType.write("  PEL   %d  \n"%(pel_count))
 # lipidType.write("  DPPC  %d  \n"%((count-lipid_count_DHPC)/2))
 # lipidType.write("  DHPC  %d  \n"%((lipid_count_DHPC)/2))
 # lipidType.write("  DPPC  %d  \n"%((count-lipid_count_DHPC)/2))
 # lipidType.write("  DHPC  %d  \n\n"%((lipid_count_DHPC)/2))
stdout.write("  Topology written in     ***pel.top*** \n")
stdout.write("  Coordinates written in  ***peg.gro***  \n\n")
