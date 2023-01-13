#! /usr/bin/env python3.6

import os 
import shutil
import numpy as np
import math
from collections import OrderedDict
import datetime
from sys import argv, stdout
from scipy import spatial

def gro_atom(line):
   """ for each line, return atom, resname, resid, x, y, z """
   # return [resid, resname, atomid, atomname, [x,y,z]]
   return(line[0:5].strip(), line[5:10].strip(), line[15:20].strip(), line[10:15].strip(),\
          [float(line[20:28]), float(line[28:36]), float(line[36:44])])

def gro_box(line):
    """ Extract the box dimension from last line of filename """
    # in general for rectangular boxes :
    # v = [float(i) for i in line.split()] + 6*[0] 
    # return v[0], v[3], v[4], v[5], v[1], v[6], v[7], v[8], v[2]

    box = [ float(i) for i in line.split()]
    return box

def read_grofile(filename):
    """ read gro file and convert all lines into list of list"""
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    return lines

def write_grofile(filename,coord, resname, atomname, box_size=None):
    """ creates a .gro file and write all the coordinates passed 
        Inputs:  filename to write coordinates
                 3x3 array of coordinates, 
                 resname,  string eg "DPPC"
                 atomname  string eg "NC3"
    """
    i = 0
    while os.path.exists(filename+"%s.gro" %i):
       i += 1       
#       j = 0 
#       while os.path.exists(filename+"%s.gro" %j):
#            j += 1
#       shutil.move(filename+"%s.gro" %i, filename+"%s.gro" %j)
 
    f = open(filename+"%s.gro" %i,'w')
    now = datetime.datetime.now().strftime("%m/%d/%Y %H:%M")
    if box_size.all() == None:
       box_size = np.array([15.0, 15.0, 15.0])
    f.write("coordinated file created by Hari on %s \n" %now)
    f.write("%d \n" %len(coord))
    for i in range(len(coord)):
        f.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' %((i+1)%100000, resname, atomname, (i+1)%100000, coord[i][0], coord[i][1], coord[i][2]))
    f.write('%8.3f%8.3f%8.3f\n' %(box_size[0],box_size[1],box_size[2]))
    f.close()
    

class Point():
      """ """
      def __init__(self, x, y, z):
          self.x = x
          self.y = y
          self.z = z

      def get_distance(self, other):
          """ Get the distance between two points"""
          dist = sqrt((self.x-other.x)**2 + (self.y -other.y)**2 +(self.z-other.z)**2)
          return dist  

      def length(self):
          length = sqrt((self.x)**2 + (self.y)**2 + (self.z)**2)
          return length           

     
      def dot_product(self, other):
          pass  

      def cross_product(self, other):
          pass
   
  
class Molecule():
    def __init__(self, filename):
# todo write docs
      self.title = ""
      self.atoms=[]
      self.box=[]
      self.filename = filename
      self.resType = ""    
      self.excludes =  []
      self.solvent  =  []
      self.coord    =  []      
      self.min_coord = []
      self.max_coord = []

      if filename:
        lines = read_grofile(filename)
        self.atoms = [gro_atom(line) for line in lines[2:-1]]
        self.box = np.array(gro_box(lines[-1]))


    def add_to_exclude(self, resType):
        """
         Add all the residues that are not needed for analysis. For eg: "W" 
         "WF"
        """ 
        self.excludes.append(resType)

    def remove_from_exclude(self, resType):
        """
         Remove the residue from exclude list that was previously.
        """
        self.excludes.remove(resType)

    def add_to_solvent(self, resType):
        """
         Add all the name of solvent (resType) in list of solvent
        """
        self.solvent.append(resType)

    def remove_from_solvent(self, resType):
        """
         Removess the residue (resType) from existing lists of solvent
        """
        self.solvent.remove(resType)

    def get_coord_all(self):
        """
            Get the coordinates of all atoms in grofile.
        """
        self.coord = [line[4] for line in self.atoms]
        return np.array(self.coord)
   
    def get_coord_all2(self, *args):
        """
            Get the coordinates of atoms with given selections.
        """

        self.coord = [line[4] for line in self.atoms]
        return np.array(self.coord)
   


    def get_box_bounds(self):
        """
           Get the min and max of box size in each x, y and z dimension
           required for dividing the box into sub-boxes.
        """
        self.get_coord_all()
        min_coord =  np.min(self.coord, axis=0)
        max_coord =  np.max(self.coord, axis=0)
        return min_coord, max_coord


    def get_coord_resType(self):
        """
            Get the coordinates of all residue types excluding solvent 
	    or any other residue type present in excludes.
        """
        selectedAtom = [line[4] for line in self.atoms if line[1] not in self.solvent and line[1] not in self.excludes]
        return np.array(selectedAtom) 

    def get_coord_resname(self, resname):
        """
           Get the coordinates of all resname eg."SOL" 
        """
        self.resname  = "resname"
        selectedAtom = [line[4] for line in self.atoms if line[1]==resname ]
        return np.array(selectedAtom)



    def get_coord_atom(self, atom, resname=None):
        """
           Get the coordinates of atom "atom" eg."NC3" 
           or the atom belonging to certain residue type (resType) eg. "DPPC".
        """
        self.atom = "atom"
        self.resname  = "resname"
        if not resname:
           selectedAtom = [line[4] for line in self.atoms if line[3]==atom]
           return np.array(selectedAtom)
        else:
           selectedAtom = [line[4] for line in self.atoms if  (line[1]==resname and  line[3]==atom)]
           return np.array(selectedAtom)


    def get_coord_atom2(self, atom, resname=None):
        """
           Get the coordinates of atom "atom" eg."NC3" 
           or the atom belonging to certain residue type (resType) eg. "DPPC".
        """
        self.atom = "atom".split()
        self.resname = "resname".split()
        print(atom)
        if not resname:
           selectedAtom = [line[4] for line in self.atoms if line[3] in atom]
           return np.array(selectedAtom)
        else:
           selectedAtom = [line[4] for line in self.atoms if  (line[1] in resname and  line[3] in atom)]
           print(self.atoms[14406])
           return np.array(selectedAtom)


    def is_keyword(self, val):
       return (val in self.kw or 
               val is None) 

    def is_keyoperation(self, val):
       return ( val in self.ko)

    kw = ['resname', 'name']
    ko = ['and', 'or']


    def select_atoms(self, selectionString, **kwargs):
        #keywords lists
        #NOT TESTED FOR DIFFERENT SITUATIONS
        """
        Only works for 4 different case example:
        select_atoms(resname ...... and/or name ....)
        select_atoms(resname)
        select_atoms(name)
        
        """
        keywords = {}
        keyoperations = {}
        keywordoperations = OrderedDict()
        keycode = {'resname': 'line[1]', 'name': 'line[3]', 'and': 'and', 'or':'or'}
        tokens = selectionString.split()

        while len(tokens) != 0:
           values = []
           if self.is_keyword(tokens[0])==True:
              word = tokens[0]
              tokens.pop(0)
              keywords[word] = []       
              keywordoperations[word]= []
           elif self.is_keyoperation(tokens[0])==True:
              operation = tokens[0] 
              tokens.pop(0)
              keyoperations[operation] = []
              keyoperations[operation] = operation
              keywordoperations[operation] = operation

          
           elif self.is_keyword(tokens[0])==False:
                val = tokens[0]
                tokens.pop(0)
                values.append(val)
                  
                keywords[word].append(val) 
                keywordoperations[word].append(val)
        a = []
        for i in range(len(list(keywordoperations.keys()))):
            a.append(list(keywordoperations.keys())[i])
     

        selectedAtom =[] 
        atomMass=[]
        mass = 'mass' in kwargs.keys()
        if mass:
           import mass_data

        if len(keywordoperations.keys())==3 and a[0]=="resname" and a[1]=="and" and a[2]=="name":
             for line in self.atoms:
                 if  (line[1] in keywords[a[0]] and line[3] in keywords[a[2]]):
                     selectedAtom.append(line[4])
#                    selectedAtom = [line[4] for line in self.atoms if  (line[1] in keywords[a[0]] and line[3] in keywords[a[2]])]
                     if mass: 
                        atType = line[3]
                        atomMass.append(mass_data.atomic_mass[atType])


        if len(keywordoperations.keys())==3 and a[0]=="resname" and a[1]=="or" and a[2]=="name":
             for line in self.atoms:
                 if  (line[1] in keywords[a[0]] or line[3] in keywords[a[2]]):
                     selectedAtom.append(line[4])
                     if mass:
                        atType = line[3]
                        atomMass.append(mass_data.atomic_mass[atType])

        if len(keywordoperations.keys())==1 and a[0]=="resname":
             for line in self.atoms:
                 if  line[1] in keywords[a[0]]:
                     selectedAtom.append(line[4])
                     if mass:
                        atType = line[3]
                        atomMass.append(mass_data.atomic_mass[atType])
 

        if len(keywordoperations.keys())==1 and a[0]=="name":
             for line in self.atoms:
                 if  line[3] in keywords[a[0]]:
                     selectedAtom.append(line[4])
                     if mass:
                        atType = line[3]
                        atomMass.append(mass_data.atomic_mass[atType])

        if mass:
           return(np.array(selectedAtom), np.array(atomMass))
        else:
           return np.array(selectedAtom)

    def get_com_all(self):
        """ 
           Get the center of mass of whole input structure.
        """
        coord = self.get_coord_all()
        com = np.mean(self.coord, axis= 0)
        return com



"""
def main():
    molecule = Molecule(argv[1])
    print(molecule.get_coord_all()) 
    print(molecule.get_box_bounds())
#    print( molecule.get_com_all())
    


if __name__=="__main__":
   main()
"""
