#! /usr/bin/env python

###############################################################
#  this script renumbers resname
#  
###############################################################

from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
from numpy import array, sum
#from termcolor import colored,cprint


resType = "EO5"
beadsResname = 9
with open('start.gro') as f1:
    with open('renumbered.gro', 'w') as f2:
      f2.write("renumbered system\n")
      f2.write("XXXXXXXX \n")
      lines = f1.readlines()
     # f2.write(lines[1])
      atomCount    = 0
      resnameCount = 0
      for i, line in enumerate(lines):
         nParticles = int(lines[1])
         j = i-2
         
         if line[5:10].strip() == resType and j%beadsResname==0:
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnameCount+1,line[5:10],line[10:15],atomCount+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
              atomCount += 1
              resnameCount+= 1
         elif i>2 and i<nParticles+2:
              f2.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %(resnameCount,line[5:10],line[10:15],atomCount+1,float(lines[i][20:28]), float(lines[i][28:36]), float(lines[i][36:44])))
              atomCount += 1
         elif i == nParticles+2:
              f2.write(line)
#              f2.write("35.000 35.000 35.000 \n")
command = " sed -i 's/XXXXXXXX/%d/g' renumbered.gro" %(atomCount)
subprocess.call(command,shell=True)
stdout.write("  Coordinates written in  ***renumbered.gro***  \n\n")
