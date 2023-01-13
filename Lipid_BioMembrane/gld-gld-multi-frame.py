#!/usr/bin/python
# written by Hari Sharma: 12/20/2015
# Calculates the center of mass of each gold and distance between gold to gold if it
# satisfies some criteria
#+++++++++++++++++++++++++++++++++++++++

from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
from numpy import array, sum
#from termcolor import colored,cprint


def read_gro(grofile):#with open('gold.gro', 'r') as grofile
   number_gold = 0
   number_of_particles = 0
   x,y,z = 0.0, 0.0, 0.0
   coordinates = []
   line_count = 0
   for line in open(grofile):
     if line_count == 1:
      number_of_particles = int(line)
      number_gold = number_of_particles/beads_gold
     elif line_count > 1 and line_count < number_of_particles + 2:
      x = float(line[20:28])
      y = float(line[28:36])
      z = float(line[36:44])
      coordinates.append([x,y,z])
     line_count += 1
   return coordinates

if len(argv) != 5:
 print '''
 ===================================================================
 Calculates the average gold to gold distance.
 Usage: You must feed the type of gold. ie c8-1nm, c16-1nm, c16-2nm.
 NOTE that index number for gold will be feed within the code. you 
 may need to change the index number for gold.
 you will need .XTC file TOPOL.TPR file

 eg: ./gld-gld.py c8-1nm traj.xtc 225000 250000
 ===================================================================\n''' 
 exit(0)
else:
 gold_ndx = 4
 gold_type = argv[1]
 trajfile = argv[2]
 initial_time = argv[3]
 final_time = argv[4]
 type = '%s' %(gold_type)
if type == 'c8-1nm':
 beads_gold = 74
 criteria = 2.88
elif type == 'c16-1nm':
 beads_gold = 110
 criteria = 4.6
elif type == 'c12-1nm':
 beads_gold = 92
 criteria =3.82

output1 = open('avg-dist.dat','w')  # to write  average gold-gold distance
output1.write('''
=======================================
 %frame number   %gold to gold distance
=======================================''')

output2 = open("without-avg-dist.dat",'w')
output2.write('''
=======================================
  %gold to gold distance
=======================================''')

# output all frame using trjconv for the last 25 ns of 
# trajectory file.
stdout.write("output all coordinate files \n")
command1 = "echo %s | trjconv_465_mpi -s topol.tpr -f %s -b %s -e %s -sep -pbc whole -o frame_dump_.gro > /dev/null" %(gold_ndx, trajfile, initial_time, final_time)
print command1
subprocess.call(command1, shell=True)

file_count = 0
while True:
  grofile = "frame_dump_" + str(file_count) + ".gro"
  if not path.isfile(grofile) or path.getsize(grofile) == 0:
     break

  #stdout.write("Taking care of snapshot %s \n " %grofile)
  
  coordinates = read_gro(grofile)
  center_of_mass = []
  #beads_gold = 74
  number_of_particles =  len(coordinates)
  data = array(coordinates)
  number_gold = number_of_particles/beads_gold
  for i in range(number_gold):
    imin = i*beads_gold
    imax = (imin + beads_gold)
    center_of_mass.append(sum(data[imin:imax],0)/beads_gold)
 
  com = array(center_of_mass)

#########################################
# Calculating distance between the center
# of mass of gold
#########################################
  import numpy as np
  distance_com = []
  #dist = array(distance_com)

  for i in range(len(com)):
    for j in range(len(com)):
       if j > i:
         distance_com.append(np.sqrt(np.sum((com[i,:] - com[j,:])**2)))         
              
# carbon-carbon bond length = 154 pm, c-atom diameter = 170 pm.
# criteria :c8-1nm = 1nm + (~ 154+170+154+170) = 1.648 nm
#           c16-1nm = 1nm +( ~ 154 + 170+154+170+154+170+154+170)= 2.296 nm
#           c16-2nm = 2nm + ( ~ 154 + 170+154+170+154+170+154+170)= 3.296 nm
#criteria = 3.3 #( twice for gold-gold distanace)

  gold_gold_distance = []
  for i in distance_com:
    if i <= criteria:
     gold_gold_distance.append(i)
     output2.write("\n %8.3f " %(i))
  avg_gold_gold_distance = sum(gold_gold_distance)/len(gold_gold_distance)
  print "avg_gold_gold_distance = %8.3f" %avg_gold_gold_distance
  output1.write(" \n    %5d            %8.3f  " %(file_count, avg_gold_gold_distance))
  #for i in gold_gold_distance:
  # output2.write("\n %8.3f " %(i))
   
#stdout.write('''\n
#++++++++++++++++++++++++++++++++++++++++++++
#Average Gold-Gold distance for %s = %8.3f
#++++++++++++++++++++++++++++++++++++++++++++\n''' %(gold_type,avg_gold_gold_distance))

  remove(grofile)
  file_count += 1
stdout.write('''
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ************************************************************
  OUTPUT WRITEN IN FILE "avg-dist.dat" & "without-avg-dist.dat"
  ************************************************************
  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n \n''')
output1.close()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

