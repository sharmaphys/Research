#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)
# Calculates the order parameter of lipids at different radius chosen at selection.dat
# It is modified from the do-order-multi.py available in martini webpage.
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



def read_gro(file,atoms):#with open('gold.gro', 'r') as grofile
   number_of_lipids = 0
   number_of_particles = 0
   x,y,z = 0.0, 0.0, 0.0
   coordinates = []
   line_count = 0
   for line in open(file):
    if line_count == 1:
     number_of_particles = int(line)
     number_of_lipids = number_of_particles/beads_lipids
    elif line_count > 1 and line_count < number_of_lipids*beads_lipids + 2:
     x = float(line[20:28])
     y = float(line[28:36])
     z = float(line[36:44])
     coordinates.append([x,y,z])
    line_count += 1

   center_of_mass = []
   data = array(coordinates)
   for i in range(number_of_lipids):
     imin = i*beads_lipids
     imax = (imin + beads_lipids)
     center_of_mass.append(sum(data[imin:imax],0)/beads_lipids)
   com = array(center_of_mass)
   
   #reading box size and making grids
   grofile = open(file) 
   grofile_lines = grofile.readlines()
   box_size = grofile_lines[-1]
   size = list(box_size.split())
#   box_x = math.ceil(float(size[0]))
#   box_y = math.ceil(float(size[1]))
   box_x = int(float(size[0])) + 1
   box_y = int(float(size[1])) + 1
#   gridWidth = 1
   grid_x = int(round(box_x/gridWidth))
   grid_y = int(round(box_y/gridWidth))
#   numberGrids = grid_x*grid_y
#   numberGrids = int(box_x/gridWidth)*int(box_y/gridWidth)
 #  numberGrids = int((box_x*box_y)/(gridWidth**2))

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
   
#   nonEmptyGrid = 0
#   for i in range(len(lipidsInGrids)):
#      if len(lipidsInGrids[i]) != 0 :
#        nonEmptyGrid += 1

           
   first   = [[] for i in range(numberGrids)]
   second  = [[] for i in range(numberGrids)]
   for i in range(numberGrids):
     if len(lipidsInGrids[i]) != 0:
       for j in lipidsInGrids[i]:
         counter = 0
         for line in open(file):
           if counter > j*beads_lipids+1 and counter < j*beads_lipids + beads_lipids+2:
              if line[10:15].strip() == atoms[0]:
               first[i].append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
              elif line[10:15].strip() == atoms[1]:
               second[i].append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
            #   print counter
           #    print second[i]
           counter += 1

   return [first, second, gridWidth, numberGrids, lipidsInGrids]

### REAL STUFF

if len(argv) != 9:
  # coments/usage
  print '''
  Compute (second rank) order parameter, defined as:

    P2 = 0.5*(3*<cos²(theta)> - 1)

  where "theta" is the angle between the bond and the bilayer normal.
  P2 = 1      perfect alignement with the bilayer normal
  P2 = -0.5   anti-alignement
  P2 = 0      random orientation

  All lipids defined in the "martini_v2.0_lipids.itp" file can be analyzed
  with this script.
  Usage: %s <traj file> <initial time> <final time> <skip time> <bilayer normal - xyz>  <lipid type>

    > %s traj.xtc 0 10000 600 0 1 0  DPPC

  will for example read a 10ns trajectory of 64 DSPC lipids, calculating the order parameter for 
  every 5th frame and averaging the results. P2 will be calculated relative to the y-axis.

  WARNING script will output all frames in one go, into files called frame_dump_XXX.gro and 
  then remove them so don't have any other files with this name in the current directory.
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
  # (normalized) orientation of bilayer normal
  orientation_of_bilayer_normal = [float(argv[5]), float(argv[6]), float(argv[7])]
  norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
  for i in range(3):
    orientation_of_bilayer_normal[i] /= norm
  stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
    orientation_of_bilayer_normal[0], \
    orientation_of_bilayer_normal[1], \
    orientation_of_bilayer_normal[2]  \
  ))
  
 
  lipid_type = argv[8]
  # output legend
  phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
  # PCs
  if   lipid_type == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  #elif lipid_type == "DHPC": bond_names = "C1B-C2B\n"
  
  elif lipid_type == "DHPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
  elif lipid_type == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
#  elif lipid_type == "DPPC": bond_names =  "C3B-C4B\n"
  elif lipid_type == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DSPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  # PEs
  elif lipid_type == "DHPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"
  elif lipid_type == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DSPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1B GL2-C1A C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  # PPCS
  elif lipid_type == "PPCS": bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 AM1-C1A GL2-D1B C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B\n"
  # output legend
  output_legend = "  Frame" + bond_names 

  # write the stuff
  stdout.write("\n " + output_legend)
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")
#  output = open('order-%s.dat' %lipid_type, 'w')
#  output.write(output_legend)
#  output.write(("-"*(len(output_legend) - 1)) + "\n")

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
  stdout.write("**********  Starting P2 calculation...  \n")
  #order_parameters = [] # for each frames
  file_count = 0
  bonds = []
  while True:
    filename = "frame_dump_" + str(file_count) + ".gro"
    if not path.isfile(filename) or path.getsize(filename) == 0:
        break
    
    stdout.write("**********  Taking care of snapshot...  %s \n\n" % filename)
    
    
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
  

    # compute order parameter for each bond, for each snapshot
    current_order_parameter = [[] for i in range(numberGrids)]
    # bonds respectively involved in the head,
    #                             in the junction head-tail,
    #                             in each tail
    bonds = []

    for bond_name in bond_names.split():
      bonds.append(bond_name.split("-"))
    for bond in bonds:
      stdout.write('*****  Taking Care of Bond : ' + str(bond) +'\n\n')
    
      # parse .gro file, grep bead coordinates
      first = read_gro(filename, bond)[0]
      second = read_gro(filename, bond)[1]
      gridWidth = read_gro(filename, bond)[2]
      lipidsInGrids = read_gro(filename, bond)[4]
      numberGrids = read_gro(filename, bond)[3]
      #stdout.write(' Lipids in  Grids: ' + str(lipidsInGrids) +'\n')
    #  stdout.write('*****  Total number of lipids in Grids: ' + str(totalLipidsInGrids) +'\n') 
      stdout.write('*****  Grids Width : ' + str(gridWidth) +'\n')
      stdout.write('*****  Number of grids : ' + str(numberGrids) +'\n\n')
      #stdout.write('  second: ' + str(second) +'\n\n\n')
     # print len(second)
      nGrids = numberGrids
      # compute order parameter for each grid
      
      order_parameter_grid = [[] for i in range(nGrids)]
      for i in range(nGrids):
        if len(first[i]) == 0:
          order_parameter_grid[i].append(0)
        elif len(first[i]) != 0:
          flist = first[i]
          slist = second[i]          
          cosThetaSquare = 0.0
          count = 0
          for j in range(len(flist)):
#            print flist[j]
             # vector between the two previous beads (orientation doesn't matter)
      #      count +=1
            vector = [0.0, 0.0, 0.0]
            for jj in range(3):
              #print flist[j]
              #print slist[j]
              vector[jj] = flist[j][jj] - slist[j][jj]
            norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
              # compute projection on the bilayer normal
            projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
              # order parameter
            cosThetaSquare += (projection**2/norm2)
          cosThetaSquare = 0.5*(3*(cosThetaSquare/len(flist))-1.0)
          order_parameter_grid[i].append(cosThetaSquare)
        current_order_parameter[i].append(order_parameter_grid[i])         

    remove(filename)
    file_count += 1
  # End while loop

  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  stdout.write("***** Snapshots analysis done.%s\n" % (" "*56))
  stdout.write("***** Computing averages...\n")
  #print current_order_parameter



  order_parameters = [[] for i in range(numberGrids)]
  for i in range(numberGrids):
   tailSn= 0.0
   eachGrid = current_order_parameter[i]
   for item in eachGrid[3:]:
     tailSn += abs(item[0])
     avgTailSn = tailSn/len(eachGrid)  
   order_parameters[i].append(avgTailSn)
   

  orderParameterMatrix = np.zeros((box_x, box_y))
#  for i in range(len(order_parameters)): 
  count = 0 
  for ii in range(0, box_x, gridWidth):
   for jj in range(0, box_y, gridWidth):
     orderParameterMatrix[ii, jj] += order_parameters[count]
     count += 1
#  print order_parameters
  print orderParameterMatrix
  np.savetxt('matrix-%s.dat'%(lipid_type), orderParameterMatrix, fmt='%8.3f',delimiter=',')
  # Write abs average order parameters <Sn> (for carbon chains only)
  # WARNING this works with currenct lipids (all have defined x5 none carbon bonds) but for manually added lipids this might not be true
  #ave_chain_s = 0
  #for i in averaged_order_parameters[3:]: 
   #  ave_chain_s += abs(i)
  #average_txt = "Abs average order parameters for carbon chains <Sn> = %8.3f \n\n" % (ave_chain_s / (len(averaged_order_parameters)-3))
  #stdout.write(average_txt)
  #output.write(average_txt)
  #stdout.write("Results written in \"grid-order-%s.dat\". \n" %lipid_type) 

  #output_order = open('order-parameter.dat', 'a')
  #output_order.write(" <Sn> of %s = %8.3f \n "%(lipid_type, ave_chain_s / (len(averaged_order_parameters)-3)))

  #output.close()
  #output_order.close()
  
  import matplotlib.pyplot as plt

  z = orderParameterMatrix
  nx,ny = np.shape(z)
  cs = plt.contourf(z)
  cb = plt.colorbar(cs, orientation = 'vertical')
  cb.set_label('P2')
  plt.xlim(0,nx)
  plt.ylim(0,ny)
  plt.savefig('contour-plot-%s.png'%lipid_type)
  plt.show()  
 

