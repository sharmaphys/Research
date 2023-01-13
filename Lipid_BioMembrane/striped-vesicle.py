#! /usr/bin/python

###############################################################
#  this script replaces the DPPC lipid by DOPC 
#  in the vesicle and generates striped vesicle
###############################################################


from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
from numpy import array, sum
#from termcolor import colored,cprint

lipidType1 = "DPPC"
beads_dppc = 12
beads_dopc = 12
lipidType2 = "DOPC"
headBead = "NC3"
stripThickness = 3

top=open('striped-vesicle.top','w')
top.write( '''
#include "martini_v2.0.itp"
#include "martini_v2.0_lipids.itp"

[ system ]
Striped-Vesicle

[ molecules ]\n
''')
with open('vesicle.gro') as f1:
    with open('striped-vesicle.gro', 'w') as f2:
      f2.write("Striped DPPPC-DOPC lipid system\n")
      
      lines = f1.readlines()
      f2.write("XXXXXXXX \n")
      total_dppc_count = 0
      total_dopc_count = 0
      totalLipidCount  = 0
      atom_count       = 0
      zCoord = []
      xCoord = []
      for line in lines:
        if line[5:10].strip() == lipidType1 and line[10:15].strip()== headBead:
          zCoord.append(float(line[36:44]))
          xCoord.append(float(line[20:28]))
      xMin  = min(xCoord)
      xMax  = max(xCoord)
      zMin  = min(zCoord)
      zMax  = max(zCoord)    
      stdout.write('z_Min= %d Z_max=%d X_Min=%d , X_Max=%d' %(zMin, zMax, xMin, xMax))
      iteration = 0
      while zMin<zMax:
        lipidCount = 0
        #dopc_count = 0
       # iteration  = 0
        for i, line in enumerate(lines):
         if line[10:15].strip() == headBead and float(line[36:44]) >= zMin and float(line[36:44]) \
          < zMin + stripThickness: 
            if iteration%2==0: 
                lipidType = lipidType1
                f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(totalLipidCount+1,lipidType,lines[i][10:15],\
                  atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
                totalLipidCount += 1
                atom_count += 1
                lipidCount += 1
              
                for j in range(1, beads_dppc):

                  f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(totalLipidCount,lipidType,lines[i+j][10:15], \
                  atom_count+1,float(lines[i+j][20:28]), float(lines[i+j][28:36]), float(lines[i+j][36:44])))
                  atom_count += 1

            elif iteration%2!=0:
                lipidType = lipidType2
                f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(totalLipidCount+1,lipidType,lines[i][10:15],\
                 atom_count+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
                totalLipidCount += 1
                atom_count += 1
                lipidCount += 1

                for j in range(1, beads_dopc):

                  f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(totalLipidCount,lipidType,lines[i+j][10:15], \
                  atom_count+1,float(lines[i+j][20:28]), float(lines[i+j][28:36]), float(lines[i+j][36:44])))
                  atom_count += 1
        top.write("%s     %d \n"%(lipidType,lipidCount))
        zMin += stripThickness
       # xMin += stripThickness
        iteration += 1    

      

      f2.write("35.000 35.000 35.000 \n")
command = " sed -i 's/XXXXXXXX/%d/g' striped-vesicle.gro" %(atom_count)
subprocess.call(command,shell=True)
#stdout.write(" \n")
#stdout.write("  DPPC  Count:-  %d  \n"%(lipid_count_DPPC))
#stdout.write("  DOPC Count:-   %d  \n\n"%(lipid_count_DOPC))
"""with open ('stripe-vesicle.top','w') as lipidType:
  lipidType.write( '''
#include "martini_v2.0.itp"
#include "martini_v2.0_lipids.itp"


 [ system ]
  Striped-Vesicle

[ molecules ]
  DPPC %d
  DPPC %d
  DHPC %d
  DPPC %d
  DHPC %d
  NA+ %d 
  ''' %(repeatUnit, pel_count, (count-lipid_count_DHPC)/2, (lipid_count_DHPC)/2, (count-lipid_count_DHPC)/2,(lipid_count_DHPC)/2,pel_count )

  )"""
stdout.write("  Topology written in     ***striped-vesicle.top*** \n")
stdout.write("  Coordinates written in  ***striped-vesicle.gro***  \n\n")
