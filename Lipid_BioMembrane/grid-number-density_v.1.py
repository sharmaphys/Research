#!/usr/bin/python
# -*- coding: utf-8 -*-

# 06/14/2016-Hari Sharma
# CALCULATES THE NUMBER DENSITY IN GRIDS using the head beads
# If you get fatal error, them take traj_dt (time) greater than how the index file are separated in radius.ndx
# for eg. if in radius.ndx, each group are separated by 300 ps time, take greater than 300 eg: dt = 400, it works.
from math import sqrt 
import math
from os import remove, system, path
from sys import argv, stdout
import subprocess
from numpy import array, sum
import numpy
import itertools
import numpy as np
numpy.set_printoptions(threshold=numpy.nan)
### SOME TOOLS
# parse a .gro file
# return a list of coordinates

def read_gro(file):
   number_of_lipids = 0
   number_of_particles = 0
   x,y,z = 0.0, 0.0, 0.0
   coordinates = []
   headBead = "NC3"
   line_count = 0
   for line in open(file):
    if line_count == 1:
     number_of_particles = int(line)
    elif line_count > 1 and line_count < number_of_particles + 2:
      if line[10:15].strip() == headBead: 
         x = float(line[20:28])
         y = float(line[28:36])
         z = float(line[36:44])
         coordinates.append([x,y,z])
         number_of_lipids += 1
    line_count += 1

   center_of_mass = []
   data = array(coordinates)
   #for i in range(number_of_lipids):
   #  imin = i*beads_lipids
   #  imax = (imin + beads_lipids)
   #  center_of_mass.append(sum(data[imin:imax],0)/beads_lipids)
   com = array(coordinates)


   import numpy as np 
   lipidsInGrids = [[] for i in range(numberGrids)]
 #  lipidsInGrids = np.empty([grid_x, grid_y])
   gridCount  = 0
   for i in range(0, box_x, gridWidth):
     for j in range(0, box_y, gridWidth):
       for lipids in range(number_of_lipids):
         if com[lipids][0] > i and  \
            com[lipids][0] < i + gridWidth and \
            com[lipids][1] > j and \
            com[lipids][1] < j + gridWidth :
             lipidsInGrids[gridCount].append(lipids)
       gridCount += 1 
   totalLipidsInGrids = sum(len(lipidsInGrids[i]) for i in range(numberGrids))

   numberDensity = [[] for i in range(numberGrids)]       
   for i in range(numberGrids):
     if len(lipidsInGrids[i]) == 0:
       numberDensity[i].append(0)
     elif len(lipidsInGrids[i]) != 0:
       numberDensity[i].append(len(lipidsInGrids[i]))

   return(numberDensity,  numberGrids, lipidsInGrids)

### REAL STUFF

if len(argv) != 6:
  # coments/usage
  print '''Usage: %s <traj file> <initial time> <final time> <skip time>   <lipid type> or both in lipid type

    > %s traj.xtc 0 10000 300 DPPC/both

  ''' % (argv[0], argv[0])
 # argv[0] is the script name in python.
  exit(0)

else:
  import numpy as np
  # snapshots
  trajfile = argv[1]
  initial_time = int(argv[2])
  final_time = int(argv[3])
  traj_dt = int(argv[4])
  lipid_type = argv[5]

  if argv[5] == both:
  	lipid_type = "DPPC DHPC"


  # Output all frame using trjconv 

  stdout.write("Output all coordinate files \n")
#  command1 = "g_select_465_mpi -sf selection.dat -selrpos whole_res_com  -f %s -s topol.tpr -on radius.ndx" % (trajfile)
  command1 = "g_select_465_mpi -select 'resname %s' -selrpos whole_res_com  -f %s -s topol.tpr -on radius.ndx" % (lipid_type, trajfile)

  print command1
  subprocess.call(command1, shell=True)
  

  m = 0
  group = 0
  for initial_time in range(initial_time, final_time, traj_dt):
   end_time = initial_time + traj_dt
#   command2 = "echo %i | trjconv_465_mpi -f %s -b %i -e %i -n radius.ndx  -dt %i  -o frame_dump_%i.gro > /dev/null" % (m, trrjfile, initial_time, end_time, traj_dt, m )
   command2 = "echo %i | trjconv_465_mpi -f %s -b %i -e %i -n radius.ndx  -dt %i  -o frame_dump_%i.gro > /dev/null" % (m, trajfile, initial_time, end_time, traj_dt, m )

   print command2
   command3 = "echo %i | editconf_465_mpi -f frame_dump_%i.gro -princ -o frame_dump_%i.gro > /dev/null" % (group, m, m)
   print command3
  # print end_time 
   command4 = "rm *# -f"
   print command4
   m += 1
   subprocess.call(command2, shell=True)
   subprocess.call(command3, shell=True)
   subprocess.call(command4, shell=True)
 

# For each dumped frame
  stdout.write("**********  Starting density calculation...  \n")
  #order_parameters = [] # for each frames
  file_count = 0
  #atoms = []
  while True:
    filename = "frame_dump_" + str(file_count) + ".gro"
    if not path.isfile(filename) or path.getsize(filename) == 0:
        break
    
    stdout.write("**********  Taking care of snapshot %s \n\n" % filename)
    
    
## reading .gro file and getting the box size, dividing into grids
    grofile = open(filename)
    grofile_lines = grofile.readlines()
    number_of_beads = grofile_lines[1]
    box_size = grofile_lines[-1]
    size = list(box_size.split())
    box_x = int(float(size[0])) + 1
    box_y = int(float(size[1])) + 1
    gridWidth = 1
#   grid_x = box_x/gridWidth
#   grid_y = box_y/gridWidth
#    numberGrids = int(box_x/gridWidth)*int(box_y/gridWidth)
    numberGrids = int((box_x*box_y)/(gridWidth**2))


    type = '%s' %(lipid_type)
    if type == 'DPPC':
     beads_lipids = 12
     number_of_lipids = int(number_of_beads)/beads_lipids
     stdout.write("*****  Lipid Type:        %s \n" % type)
     stdout.write("*****  Number of Lipids:  %s \n\n" % number_of_lipids)
    elif type  == 'DHPC':
     beads_lipids = 8
     number_of_lipids = int(number_of_beads)/beads_lipids
     stdout.write("*****  Lipid Type:        %s \n" % type)
     stdout.write("*****  Number of Lipids:  %s \n\n" % number_of_lipids)
    elif type == 'DPPC DHPC':
     stdout.write("***** Calculated density of Lipid Type:        %s \n" % type)

#    stdout.write('*****  Taking Care of atom : ' + str(atom) +'\n\n')
    
      # parse .gro file, grep bead coordinates
    (numberDensity, lipidsInGrids,numberGrids) = read_gro(filename)
    #lipidsInGrids = read_gro(filename)[2]
    #numberGrids = read_gro(filename)[1]
    stdout.write('*****  Grids Width : ' + str(gridWidth) +'\n')
    stdout.write('*****  Number of grids : ' + str(numberGrids) +'\n\n')
    nGrids = numberGrids

    remove(filename)
    file_count += 1
  # End while loop

  #stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  stdout.write("***** Snapshots analysis done.%s\n" % (" "*56))
  stdout.write("***** Computing averages...\n")

  densityMatrix = np.zeros((box_x, box_y))
  count = 0 
  for ii in range(0, box_x, gridWidth):
   for jj in range(0, box_y, gridWidth):
     densityMatrix[ii, jj] += numberDensity[count]
     count += 1
  #print densityMatrix
  np.savetxt('density-matrix-%s.dat'%(lipid_type), densityMatrix, fmt='%8.3f',delimiter=',')
  
  import matplotlib.pyplot as plt

  z = densityMatrix
  nx,ny = np.shape(z)
  cs = plt.contourf(z)
  cb = plt.colorbar(cs, orientation = 'vertical')
  cb.set_label('number density')
  plt.xlim(0,nx)
  plt.ylim(0,ny)
  plt.savefig('density-plot-%s.png'%lipid_type)
  plt.show()  
 

