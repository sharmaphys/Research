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
     elif line_count > 1 and line_count < number_gold*core + 2:
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

if type == 'c4-1nm':
 beads_gold = 56 
 criteria = 1.94
 core = 38
if type == 'c8-1nm':
 beads_gold = 74
 criteria = 2.88
 core = 38
if type == 'c12-1nm':
 beads_gold = 92
 criteria = 3.82
 core = 38
if type == 'c16-1nm':
 beads_gold = 110
 core = 38
 criteria = 4.76
elif type == 'c16-2nm':
 beads_gold = 212 # only gold core
 core = 140
 criteria = 6.6

 # carbon-carbon bond length = 0.47 nm in MARTINI  ##############################
# criteria :c8-1nm = 0.5 + (0.47*2)  = 1.44 nm >> criteria = 2*1.44 = 2.88 nm  #
# criteria :c12-1nm = 0.5 + (0.47*3) = 1.91 nm >> criteria = 2*1.91 = 3.82 nm  #
# criteria :c16-1nm = 0.5 + (0.47*4) = 2.38 nm >> criteria = 2*2.38 = 4.76 nm  #
################################################################################

output1 = open('avg-core-dist.dat','w')  # to write  average gold-gold distance
output1.write('''
===============================================
 %frame number % avgerage core to core distance
===============================================''')

output2 = open('core-core-dist.dat','w')
output2.write('''
=======================================
  % core to core distance
=======================================''')

output3 = open('cluster-size.dat','w')  # to write number of gold clusters
output3.write('''
===============  =======================  =================
  frame_number   maximum_cluster_size     number_of_clusters
===============  =======================  =================''')

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
  command2 = "sed -i '/CCC/d' %s " %(grofile)
#  print command2
  subprocess.call(command2, shell=True)
  command3 = "sed -i '/SSS/d' %s " %(grofile)
  #print command3
  subprocess.call(command3, shell=True)
  coordinates = read_gro(grofile)
  center_of_mass = []
  #beads_gold = 74
  number_of_particles =  len(coordinates)
  data = array(coordinates)
  number_gold = number_of_particles/core
  for i in range(number_gold):
    imin = i*core
    imax = (imin + core)
    center_of_mass.append(sum(data[imin:imax],0)/core)
 
  com = array(center_of_mass)

#########################################
# Calculating distance between the center
# of mass of gold
#########################################
  import numpy as np
  cluster_matrix = np.zeros((number_gold,number_gold)) ## to get matrix with clustered gold
  distance_com = []
  gold_gold_distance = [] 
  for i in range(len(com)):
    for j in range(len(com)):
       if j > i :  
         dist=(np.sqrt(np.sum((com[i,:] - com[j,:])**2)))
         if dist <= criteria:
            
  #gold_gold_distance = []
  #for i in distance_com:
  #  if i <= criteria:
          gold_gold_distance.append(dist)

  #for i in range(len(gold_gold_distance)):
  
          cluster_matrix[i][j] += 1  # to get the matrix of clustered gold nanoparticles
          output2.write("\n %8.3f " %(dist))

  avg_gold_gold_distance = sum(gold_gold_distance)/len(gold_gold_distance)
  print "\n  avg_gold_gold_distance = %8.3f" %avg_gold_gold_distance
  output1.write(" \n        %5d          %8.3f  " %(file_count, avg_gold_gold_distance))
  #print cluster_matrix


  connected_gold_id = []

  for i in range(number_gold):
   connected_gold_id.append([])
   for j in range(number_gold):
     if  cluster_matrix[i][j] == 1:
       for k in range(number_gold): 
        if i not in connected_gold_id[i]:
         connected_gold_id[i].append(i)
        if j not in connected_gold_id[i]:
         connected_gold_id[i].append(j)   
  #print connected_gold_id
  
  cluster = connected_gold_id
  l = len(cluster)
  i = 0
  while i < (l-1):
      for j in range(i+1, l):

       # i,j iterate over all pairs of l's elements including new 
       # elements from merged pairs. We use len_l because len(l)
       # may change as we iterate 
       i_set = set(cluster[i])
       j_set = set(cluster[j])
  
       if len(i_set.intersection(j_set)) > 0:
      
         # Remove these two from list
         cluster.pop(j) 
         cluster.pop(i)

         # Merge them and append to the orig. list
         ij_union = list(i_set.union(j_set))
         cluster.append(ij_union)

         # len(cluster) has changed
   
         l -= 1

         # adjust 'i' because elements shifted
         i -= 1

         # abort inner loop, continue with next l[i]
         break
   
      i += 1
  cluster.sort()
 ## remove empty lists
  cluster = [x for x in cluster if x]
#  print cluster
  
  cluster_size = []
  for i in range(len(cluster)):
    cluster_size.append(len(cluster[i]))
  #print cluster_size
  not_clustered_gold = number_gold - sum(cluster_size)
  #print "\n  Not clustered gold = %i" %not_clustered_gold
  
  max_cluster_size = max(cluster_size)
  number_cluster   = len(cluster_size) + not_clustered_gold
  output3.write("\n   %5d              %5d                   %5d" %(file_count, max_cluster_size, number_cluster))
              
  remove(grofile)
  file_count += 1
stdout.write('''
  =======================================================================================
  OUTPUT WRITEN IN FILE:: "avg-gold-gold-dist.dat", "core-core-dist.xvg" & "cluster-size.dat"
  =======================================================================================\n \n''')
output1.close()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
