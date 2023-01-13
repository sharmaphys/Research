#! /usr/bin/env python
"""
   Description:
   --------------------------------------------------------------------------------------------------
    Calculates the persistence length of the polymer chain.
    You will need traj.xtc (or .gro)  and topol.tpr file.
   
    Usage: ./script_name -f traj.xtc -n index_number -s topol.tpr -b 0 -e 500
           ./script_name -g filename.gro

    Example 1: to calculate the persitence length of a polymer chain 
              having resname ST1 ST2 ST3 with backbone atoms C11 C21 C31, use

             ./script_name -g filename.gro -lp "resname ST1 ST2 ST3 and name C11 C21 C31"

    Example 2: to calculate the radius of gyration of a polymer chain 
              having resname ST1 ST2 ST3 , use

             ./script_name -g filename.gro -rg "resname ST1 ST2 ST3"

    Example 3: to calculate the end to end distance of a polymer chain 
              having resname ST1 ST2 ST3 with ST1 and ST3 terminal groups and 
              atoms C31 in ST1 and C31 in ST3 as terminal atoms, use

             ./script_name -g filename.gro -ee "resname ST1  and name C31" "resname ST3  and name C31"

    Example 4: to calculate all Rg, end to end distance and persistence length of a polymer chain 
              having resname ST1 ST2 ST3 and backbone atoms C11 C21 and C31, use

             ./script_name -g filename.gro  -lp resname "ST1 ST2 ST3 and name C11 C21 C31"
             -rg "resname ST1 ST2 ST3" -ee "resname ST1  and name C31" "resname ST3  and name C31"
     
    Example 5: USING Trajectory file (.xtc) to calculate all Rg, end to end distance and persistence 
               length of a polymer chain having resname ST1 ST2 ST3 and backbone atoms C11 C21 and C31, use

             ./script_name -f traj.xtc -s topol.tpr -n index number -b (beginning time) -e (end time) 
             -lp "resname ST1 ST2 ST3 and name C11 C21 C31"
             -rg "resname ST1 ST2 ST3" -ee "resname ST1  and name C31" "resname ST3 and name C31"

    Example 5: USING Trajectory file (.xtc) to calculate segmental autocorrelation, use

             ./script_name -f traj.xtc -s topol.tpr -n index number -b (beginning time) -e (end time) 
             -ee "resname ST1 and name C31 " "resname ST3 and name C31 " 
             -segacf "resname ST1 ST2 ST3 and name C31 C21 C11" 5
           
             # note that currerntly at least one of the -ee or -lp or -rg are required.
             # need to make it independent
   ---------------------------------------------------------------------------------------------------------
"""
# Written by Hari Sharma: 08/06/2018
# version 0.2 

from math import sqrt, log
from os import remove, system, path
from sys import argv, stdout
import subprocess
from statistics import stdev
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy
from numpy import array, sum
import numpy as np
from numpy import linalg as LA
import os
import timeit
from itertools import chain
import structure, topology, mass_data

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()

def mean(x):
   """"
    Returns the arithematic mean of values in 1-d array
   """
   n = len(x)
   if n < 1:
      raise ValueError('NOt sufficient data. Mean requires at least one data') 
   avg = sum(x)/float(n)
   return avg

def standard_dev(x):
   """
   Get the standard deviation of the 1-d array.
   """
   n = len(x)
   if n <2:
      raise ValueError(' variance requires at least two data points')
   c = mean(x)  
   variance = sum((x -c)**2)/float(n)
   sd = sqrt(variance)
   return sd 

def read_gro(grofile):
   """
   Reads the coordinate(.gro) file and get the coordinates
   """
   nAtoms = 0
   x,y,z = 0.0, 0.0, 0.0
   coordinates = []
   lineCount = 0
   for line in open(grofile):
     if lineCount == 1:
      nAtoms = int(line)
     elif lineCount > 1 and lineCount < nAtoms + 2: 
              x = float(line[20:28])
              y = float(line[28:36])
              z = float(line[36:44])    
              coordinates.append([x,y,z])
     lineCount += 1
   coordinates = array(coordinates)
   N = len(coordinates)
   return(coordinates, N)


def get_bond_vec(coordinates):
    """
    Get the normalized vector between two points and its norm
    given a array of coordinates.
    """
    N      = len(coordinates)
    nBonds = N-1
    bondVector = []
    norm       = np.zeros(nBonds)
    x, y, z = 0.0, 0.0 , 0.0
    for i in range(nBonds):
         x = coordinates[i+1][0] - coordinates[i][0]
         y = coordinates[i+1][1] - coordinates[i][1]
         z = coordinates[i+1][2] - coordinates[i][2]
         norm[i] += sqrt(x**2 + y**2 + z**2)
         bondVector.append([x/float(norm[i]),y/float(norm[i]),z/float(norm[i])])
    bondVector = array(bondVector)
    return (bondVector, norm)


def traj_file():
    """
    Extract .gro files from a .xtc file from initial
    time to final time and do rest of the analysis.

    """
    if not arg.index_no:
       print('\t \n What is the index number for the group you want? \n')
       exit(0)
    if not arg.topol_tpr:
       print('\t \n .tpr file needed. \n')
       exit(0)
    if not arg.begin_time:
       print('\t \n Initial time= ?? \n')
       exit(0)
    if not arg.end_time:
       print('\t \n Final time= ?? \n')
       exit(0)
 
    xtc = arg.traj_file
    tpr = arg.topol_tpr
    t1  = arg.begin_time
    t2  = arg.end_time
    n   = arg.index_no 

    command1 = "echo %s | trjconv_mpi -s %s -f %s -b %s -e %s -n -sep -pbc whole -o frame_dump_.gro > /dev/null"  %(n, tpr, xtc, t1, t2)
    subprocess.call(command1, shell=True)
    return None


def persistence_ete_projection(coordinates):
    """
    Calculates the local persistence length along the backbone 
    using the projection of the backbone end-to-end vector
    Re on the generic kth spring.(doi: 10.1063/1.1651052)

    lp(k) = < (rk,k+1/(|rk,k+1|)).Re >

    Also, calculates the persistence length using Flory's 
    approach using the projections of all bond vectors on the 
    first vector and is considered to be the exact and without 
    any approximation (10.1016/j.polymer.2004.06.034).

    In general, both approaches should be similar, however turns 
    out to be different in case of not having approximately equal
    bond lengths. (don't know why) 

    Exact persistence length is Lp = < sum over bonds ri.r1 >/ <r>, 
    where <r> is the average bond length, r1 is the first bond 
    vector and ri are all bond vectors.ie average sum of the 
    projection of backbone end to end vector on the first vector

    """
    nBonds = len(coordinates)-1 # no of bonds
    bondVector = []
    bond_vector_normalized = []
    norm       = np.zeros(nBonds)
    x, y, z = 0.0, 0.0 , 0.0
    for i in range(nBonds):
         x = coordinates[i+1][0] - coordinates[i][0]
         y = coordinates[i+1][1] - coordinates[i][1]
         z = coordinates[i+1][2] - coordinates[i][2]
         norm[i] += sqrt(x**2 + y**2 + z**2)
         bond_vector_normalized.append([x/float(norm[i]),y/float(norm[i]),z/float(norm[i])])
         bondVector.append([x, y, z])

    coordinates = np.array(coordinates)
    # end to end vector
    re = coordinates[-1]-coordinates[0]

    # average bond length
    avg_bond_length = np.mean(norm)

    # for projection of Ree on individual bond
    lp_k = []
    lp_exact = 0
    for i in range(nBonds):
        lp_k.append(np.dot(bond_vector_normalized[i], re))
#        lp_k.append(np.inner(bondVector[i], re)/avg_bond_length)

        # note that we are dividing by average bond length, will get same result if use norm[0]
        lp_exact += float(np.inner(bondVector[0], bondVector[i])/avg_bond_length)
       
    return np.absolute(lp_k), np.absolute(lp_exact)


def get_segment_bond_vector(coordinates, sep):
    """
    Given the coordinates of (backbone), get the vector
    connecting between atoms separated by sep bonds.
    """
    N = len(coordinates)
    if sep > N:
       print("\t Bond vector separation greater than the number of atoms in backbone. \n")
       exit(0)
    nSegments = int(N/sep)
    bondVector = []
    for i in range(nSegments):
        ii = i * sep
        jj = (ii + sep) -1  # -1 is to get the correct indicies
        vec = coordinates[jj]-coordinates[ii]
        bondVector.append([vec[0], vec[1], vec[2]])
    return bondVector   # return as list 

    

def vector_autocorrelation(eteVector):
    """
    Calculates the normalized second rank end-to-end vector autocorrelation.
     P1 = <cos(theta)>
     P2 = 0.5*(3<cos(theta)> -1) , where cos(theta) is the angle between the 
     end to end unit vectors at time t= 0 and at times t = t. 

    input: of end to end vector, nx3 

    returns :  array of normalized autocorrelation 
    """
    N = len(eteVector)
    acf_first = []  # p1
    acf_second = [] # p2
    for i in range(N):
        dotProd = np.dot(eteVector[0], eteVector[i])
        norm0 = np.linalg.norm(eteVector[0])
        normi = np.linalg.norm(eteVector[i])
       
        cos_theta = dotProd/float(norm0*normi)
        p2 = 0.5 *( 3*(cos_theta)**2 -1)
        acf_first.append(p2)
#        acf_second.append(p2)
       
#    acf_temp = [[] for i in range(N)]
#    for i in range(N):
#        k = 0
#        for j in range(i, N):
#            dotProd = np.dot(eteVector[i], eteVector[j])
#            norm0 = np.linalg.norm(eteVector[i])
#            normi = np.linalg.norm(eteVector[j])
#            acf_temp[k].append(dotProd/(norm0*normi))
#            k += 1
#    acf = [[] for i in range(N)]
#    for i in range(N):
#        acf[i] = sum(acf_temp[i])/float(len(acf_temp[i]))

    return(np.array(acf_first))

def correlation_function_with_first(coordinates):
    """
     Calculate the Cos(theta) between the normalized vectors.
     Specifically, dot product of first vector with all the rest
     of the vectors.

     <u(0).u(i)> = <cos(theta)> 
    """
    N = len(coordinates)
    (bondVector, norm)    = get_bond_vec(coordinates)
#    coordinates = np.array(coordinates) 
    nBonds = N-1

    # for correlation only with first segment.
    CL = []
    CF = []
    s = 0
    k = 0
    cFunc = 0
    n = 1
    for i in range(0, nBonds, 1):
        s += norm[i]
#        cFunc = np.inner(bondVector[0], bondVector[i])
        cFunc += np.inner(bondVector[0], bondVector[i])
        CL.append(s)
        CF.append(float(cFunc/n))
        n += 1
    return(CL, CF)

def bond_angle_correlation(coordinates):
    """
    calculates the average of cos(theta) between two consecutive bond
    vectors along the chain;
    Then the persiscence length lp = <l>/(1-<cos(theta)>)

    """
    N = len(coordinates)
    (bondVector, norm)  = get_bond_vec(coordinates)
    nBonds = N-1
   
    cos_theta = 0.0
   
    avg_bond_length = np.mean(norm, axis= 0)    
    count = 0
    for i in range(nBonds-1):
        cos_theta += np.inner(bondVector[i], bondVector[i+1])
        count += 1
    avg_cos_theta = float(cos_theta/count) 
    lp = float(avg_bond_length/(1-avg_cos_theta))
    print("avg_bond_length", avg_bond_length, "<cos>", avg_cos_theta, "lp", lp )
    return lp

def correlation_function(coordinates):
    """
     Calculate the Cos(theta) between the normalized vectors.
     Specifically, dot product of first vector with all the rest
     of the vectors, dot product of second vector and all the
     rest of vectors except first vector and so on.
    """
    N = len(coordinates)
    (bondVector, norm)    = get_bond_vec(coordinates)    
    nBonds = N-1
    contourLength = [[] for i in range(nBonds)]
    correlationFunction = [[] for i in range(nBonds)]
    
    """
    # for correlation only with first segment.
    CL = []
    CF = []    
    s = 0
    k = 0
    cFunc = 0
    for i in range(nBonds):
        s += norm[i]
        cFunc = np.inner(bondVector[0], bondVector[i])
        CL.append(s)
        CF.append(cFunc)

    """ 
    for i in range(nBonds): 
        k = 0 
        s = 0                      
        cFunc = 0
        n = 1
        # taking steps of 2 for alternate bonds, as only trans bonds alligns(similar to 
        # how gromacs calculate.
        for j in range(i, nBonds, 1):
           # get the contour length(total bond length from bond i to bond j included, (not distance from i to j))
           s = sum(norm[i:j+1])
           cFunc += np.inner(bondVector[i], bondVector[j])
#           cFunc = np.inner(bondVector[i], bondVector[j])

           cFunc_avg = float(cFunc/n)
           contourLength[k].append(s)
           correlationFunction[k].append(cFunc_avg)
#           correlationFunction[k].append(cFunc)

           k += 1
           n += 1
    
    # remove emtpy lists 
    CL = [x for x in contourLength if x]
    CF = [x for x in correlationFunction if x]
  

    nn = len(CL)
    if len(CL)!= len(CF):
       raise ValueError('Different Length of Correlation function and Contour length')

    avgCL = [[] for i in range(nn)]
    avgCF = [[] for i in range(nn)]
    for i in range(nn):
      avgCL[i] = sum(CL[i])/float(len(CL[i]))
      avgCF[i] = sum(CF[i])/float(len(CF[i]))

#    # only positive values of correlation function, assuming it should decay to 0
#    avgCLPos = []
#    avgCFPos = []
    
#    for i in range(len(avgCL)):
#        if avgCF[i]>0:
#              avgCLPos.append(avgCL[i])
#              avgCFPos.append(avgCF[i])
#        else:
#           break
#    print("n coordinates", N, "nbonds", nBonds, len(avgCLPos))
#    return(avgCLPos, avgCFPos)

    return(avgCL, avgCF)

def interpolate(x,y):
    """
    Interpolate linear on a log scale on y.
    Given y , find x. 
    Input:
    -----   lists, x and y values
   
    Returns: Interpolated value, float  

    """
    j = 0
    for i in range(len(y)):
        if y[i] < 0.367:  # 1/e
           j += i
           break
    if j == 0:
       xx = x[i]  # persistence length for straight, rigid case is the 
                  # whole contour length

    y2 = log(y[j])
    y1 = log(y[j-1])
    yy = -1
  
    x2 = x[j]
    x1 = x[j-1]

    xx = x[j-1] + (x[j]-x[j-1])*(( 1 + log(y[j-1]))/(float(log(y[j-1]) -log(y[j]))))
        
    return xx


def fit_exponential_decay(x, y):
    r"""Fit a function to an exponential decay
    .. math::  y = \exp(-x/a)

    Parameters
    ----------
    x, y : array_like
      The two arrays of data
    Returns
    -------
    a : float
      The coefficient *a* for this decay
    Notes
    -----
    This function assumes that data starts at 1.0 and decays to 0.0
    Requires scipy
    """
    from scipy.optimize import curve_fit

    # get only positive values
    xPos = []
    yPos = []
 
    for i in range(len(x)):
        if y[i]>0:
              xPos.append(x[i])
              yPos.append(y[i])
        else:
           break

    xPos = np.array(xPos)
    yPos = np.array(yPos)
    def expfunc(x, a):
        return np.exp(-x/a)

    a = curve_fit(expfunc, xPos, yPos)[0][0]
    return a

def fit_linear(x,y):
    from scipy.optimize import curve_fit

def avgUnequalArray(unequalArray):
    """
    Returns nx2 array of average of columns(i.e, along axis=0) and its standard deviation 
    from an array of arrays having unequal length.
    Example: if a = [[1],
                     [1,2,3],
                     [2,3,4,9]], average of array(a) returns:
                     [1.33333333 2.5  3.5   9.  ]
    """
    output = []
    maxlength = 0
    for arr in unequalArray:
        if len(arr)>maxlength:
           maxlength = len(arr)

    for index in range(maxlength):
        temp = []
        for arr in unequalArray:
            if index < len(arr):
               temp.append(arr[index])
        if len(temp)>1:
           stdev = standard_dev(np.array(temp))
        else:
           stdev = 0
        output.append([mean(temp), stdev])
    return np.array(output)


def main():
    if len(argv) < 2:
       print(parser.print_help())
       exit(1)

    if arg.gro_file:
       nFrame = 1
       grofile = arg.gro_file

       output = open('cf_vs_cl_polymer.xvg','w')
       output.write("@    title \"Correlation Function\" \n")
       output.write("@    xaxis  label \"Contour Length (s)(nm)\" \n")
       output.write("@    yaxis  label \" C = <u(0).u(s)> \"\n")
       output.write("@TYPE xy""\n")

 
       results = {}
       rg = 0.0
       anisotropy = 0.0
       asphericity = 0.0
       ete = 0.0
       lp = 0.0
       lp_k = []
       lp_exact = 0.0

       if arg.rg:
          molecule = structure.Molecule(filename=grofile)
          selection_string1 = arg.rg
          (coord, mass) = molecule.select_atoms(selection_string1,mass=True)
#          rg += topology.radius_of_gyration(coord, mw = mass)
          rg, anisotropy, asphericity = topology.anisotropy(coord, mw = mass)
          results['rg'] = rg
          results['anisotropy'] = anisotropy
          results['asphericity']= asphericity


       if arg.ee:
          molecule = structure.Molecule(filename=grofile)
          selection_string2 = arg.ee[0]
          selection_string3 = arg.ee[1]
          start = molecule.select_atoms(selection_string2)
          end = molecule.select_atoms(selection_string3)
          ete += topology.end_to_end_distance(start, end)
          results['EtE'] = ete
          
       if arg.lpk:
          molecule = structure.Molecule(filename=grofile)
          selection_string4 = arg.lpk
          coordinates = molecule.select_atoms(selection_string4)
#          lp_k.append(persistence_ete_projection(coordinates)[0])
          lp1, lp2 = persistence_ete_projection(coordinates) 
          lp_k.append(lp1)         
          lp_exact += lp2
          lp_k = np.array(lp_k)
          ete_correlation = open('end_to_end_correlation_polymer.xvg', 'w')
          ete_correlation.write("@    title \"Persistence Length lp(k) vs k\" \n")
          ete_correlation.write("@    xaxis  label \"k \" \n")
          ete_correlation.write("@    yaxis  label \" Persistence Length lp(k) \"\n")
          ete_correlation.write("@TYPE xy""\n")
          for i in range(len(lp_k[0])):
              ete_correlation.write(" %8.3f    %8.3f\n" %(i+1, lp_k[0][i]))

          results['Lp_Flory'] = lp_exact         


       if arg.lp:
          molecule = structure.Molecule(filename=grofile)
          selection_string5 = arg.lp
          coordinates = molecule.select_atoms(selection_string5)
          (CL, CF) =  correlation_function_with_first(coordinates)
#          (CL, CF) =  correlation_function(coordinates)


          for i in range(len(CL)):
              output.write(" %8.3f    %8.3f\n" %(CL[i], CF[i]))

          lp += fit_exponential_decay(CL, CF)
          results['Lp_Exp'] = lp

       keys = list(results.keys())
       nKeys = len(list(results.keys()))
       keyLength = [(len(key)+3)*'=' for key in keys]
       output1 = open('lp-polymer.dat', 'a')
       output1.write(('    {}' + '\t\t{:>s}(nm)'*nKeys +'\n').format("Frame", *keys))
       output1.write(('    {:=^5}'+'\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       values = []
       for i in range(nKeys):
           values.append(results[keys[i]])
       output1.write(('    {:^5}'+'\t\t{:.3f}'*nKeys +'\n').format(i+1, *values))

       k = results.keys()
       v = results.values()
       merged_all = [j for i in zip(k,v) for j in i]
       print (("\n \t {:5s} : {:.3f} (nm)"*nKeys +'\n').format(*merged_all))
       print("""
       \t Note that only +ve values of averages of cosines are used to calculate persistence 
       \t length using exponential fit. \n""")


##       print ("\n \t Rough estimate of persistence Length using interpolation = %5.3f nm\n" %lp_interpolation)

    elif arg.traj_file:
       if arg.traj_file=='multi':
         pass
       else:# arg.traj_file == 'traj'

          traj_file() 
       nFrame = 0

       results = {} 
       
       if arg.rg:
            results['rg'] = []
            results['anisotropy'] = []
            results['asphericity'] = []
       if arg.ee:
            results['ete'] = []          
       if arg.lp:
            results['Lp_Exp'] = []
       if arg.lpk:
            lpk_total = []
            results['Lp_Flory'] = []
           
       contourLenTotal = []
       corrFuncTotal = []
       lpTotalAllFrame = 0
       eteVecTot = []
       # bond vector of segments from all frames
       segmentVector = []

       while True:
         grofile = "frame_dump_" + str(nFrame) + ".gro"
         if not path.isfile(grofile) or path.getsize(grofile) == 0:
            break


         if arg.rg:
            molecule = structure.Molecule(filename=grofile)
            selection_string1 = arg.rg
            (coord, mass) = molecule.select_atoms(selection_string1,mass=True)
#            rg = topology.radius_of_gyration(coord, mw = mass)
            rg, anisotropy, asphericity = topology.anisotropy(coord, mw = mass) 
            results['rg'].append(rg)
            results['anisotropy'].append(anisotropy)
            results['asphericity'].append(asphericity)
            

         if arg.ee:
            molecule = structure.Molecule(filename=grofile)
            selection_string2 = arg.ee[0]
            selection_string3 = arg.ee[1]
            start = molecule.select_atoms(selection_string2)
            end = molecule.select_atoms(selection_string3)              
            ete = topology.end_to_end_distance(start, end)
            results['ete'].append(ete)

         if arg.eeacf:
            molecule = structure.Molecule(filename=grofile)
            selection_string2 = arg.eeacf[0]
            selection_string3 = arg.eeacf[1]
            start = molecule.select_atoms(selection_string2)
            end = molecule.select_atoms(selection_string3)
            eteVec = (end-start)
            eteVecTot.append(list(eteVec[0])) # 0 index because it has only one row

         if arg.segacf:
            molecule = structure.Molecule(filename=grofile)
            selection_string4 = arg.segacf[0]
            separation = int(arg.segacf[1])
            atoms = molecule.select_atoms(selection_string4)
            bondVec = get_segment_bond_vector(atoms, separation) 
            segmentVector.append(bondVec)


         if arg.lp:
            molecule = structure.Molecule(filename=grofile)
            selection_string6 = arg.lp            
            coordinates = molecule.select_atoms(selection_string6)
#            (CL, CF) =  correlation_function_with_first(coordinates)
            (CL, CF) =  correlation_function(coordinates)


            contourLenTotal.append(CL)
            corrFuncTotal.append(CF)


#            for i in range(len(CL)):
#               output.write(" %8.3f    %8.3f\n" %(CL[i], CF[i]))

            lp = fit_exponential_decay(CL, CF)
            print('Frame {}, Persistence Length: {:8.3f} nm'.format(nFrame, lp))
            results['Lp_Exp'].append(lp)

         if arg.lpk:
            molecule = structure.Molecule(filename=grofile)
            selection_string7 = arg.lpk
            coordinates = molecule.select_atoms(selection_string7)
            lp_k, lp_1 = persistence_ete_projection(coordinates)     
            lpk_total.append(lp_k)
            results['Lp_Flory'].append(lp_1)

#         remove(grofile)
         nFrame += 1

       if arg.eeacf:
          acf = vector_autocorrelation(eteVecTot)
          ete_acf = open('acf.xvg', 'w')
          ete_acf.write("@    title \"End to End Vector auto correlation\" \n")
          ete_acf.write("@    xaxis  label \"time(t) \" \n")
          ete_acf.write("@    yaxis  label \" r(t).r(0) \"\n")
          ete_acf.write("@TYPE xy""\n")
          for i in range(len(acf)):
              ete_acf.write(" %5d    %8.3f \n" %(i, acf[i]))
          ete_acf.close()
 
       if arg.segacf:
          tRestart = arg.tres
          frameSep = arg.dt

          # number of frame after which to restart
          frameRestart = int(tRestart/frameSep)
          len_first_row = len(segmentVector[0])

#         acf_segments_frame = np.zeros((len(segmentVector), len(segmentVector[0]))) 
          acf_total = [] # from all frames with multiple restart points

          if frameRestart >= len(segmentVector):
             ii = 1 # to iterate over frames, 1 gives only one restart 
             jj = len(segmentVector) # to save p2 values per restart
             frameRestart = len(segmentVector)
             print(ii, jj)
          else: 
             ii = len(segmentVector) - frameRestart
             jj = frameRestart


          for i in range(ii):
#              jj = frameRestart
#              ii = len(segmentVector) -i 
              acf_segments_frame = np.zeros((jj, len(segmentVector[0])))
              icount = 0
              for j in range(i, i+frameRestart):
#              for j in range(i, len(segmentVector)):
                  for k in range(len_first_row):
                     if len(segmentVector[j][k]) !=0:
                        dotProd = np.dot(segmentVector[i][k], segmentVector[j][k])
                        norm0 =   np.linalg.norm(segmentVector[i][k])
                        norm1 =   np.linalg.norm(segmentVector[j][k])
                        cos_theta = dotProd/float(norm0*norm1)
                        p2 = 0.5*(3*(cos_theta)**2 -1)
                        acf_segments_frame[icount][k] +=p2
                  icount += 1

              acf_segment_avg = np.mean(acf_segments_frame, axis = 1)
              acf_total.append(acf_segment_avg)
#          acf_segment_avg = np.mean(acf_segments_frame[:, 1:-1], axis = 1)

          acf_total_avg = avgUnequalArray(acf_total)

          seg_acf = open('segacf.xvg', 'w')
          seg_acf.write("@    title \"Segmental auto correlation\" \n")
          seg_acf.write("@    xaxis  label \"time(t) \" \n")
          seg_acf.write("@    yaxis  label \" P2 \"\n")
          seg_acf.write("@TYPE xy""\n")
          for i in range(len(acf_total_avg)):
              seg_acf.write(" %5d    %8.3f \n" %(i, acf_total_avg[i][0]))
#              seg_acf.write((" {:>3d}" + "{:>8.3f}" +   "{}" +"\n").format(i, acf_total_avg[i], " ".join(format(x, "8.3f") for x in acf_segments_frame[i])))
          seg_acf.close() 

       if arg.lpk:
          lpk_total = np.array(lpk_total)
          ete_correlation = open('end_to_end_correlation_polymer.xvg', 'w')
          ete_correlation.write("@    title \"Persistence Length lp(k) vs k\" \n")
          ete_correlation.write("@    xaxis  label \"k \" \n")
          ete_correlation.write("@    yaxis  label \" Persistence Length lp(k) \"\n")
          ete_correlation.write("@TYPE xy""\n")
          ete_correlation.write("@ view 0.15, 0.15, 0.75, 0.85 ""\n")
          ete_correlation.write("@ legend on ""\n")
          ete_correlation.write("@ legend box on ""\n")
          ete_correlation.write("@ legend loctype view ""\n")
          ete_correlation.write("@ legend 0.78, 0.8 ""\n")
          ete_correlation.write("@ legend length 2 ""\n")
          ete_correlation.write("@ s0 legend \"Lp(k)\" ""\n")
          ete_correlation.write("@ s0 legend \"std-dev\" ""\n")



#          lpk_avg = np.mean(lpk_total, axis = 0)
#          lpk_std = np.std(lpk_total, axis = 0)

          lpk_avg = avgUnequalArray(lpk_total)
          for i in range(len(lpk_avg)):
#              ete_correlation.write(" %8.3f    %8.3f    %8.3f \n" %(i+1, lpk_avg[i] , lpk_std[i]))
              ete_correlation.write(" %8.3f    %8.3f   %8.3f \n" %(i+1, lpk_avg[i][0] , lpk_avg[i][1]))

          ete_correlation.close()


         # get the average of contour length and correlation function from all frames
       if arg.lp:
         contourLen_tot = np.array(contourLenTotal)
         corrFunc_tot = np.array(corrFuncTotal)

         cfVsCl = open('cf_vs_cl_polymer_total.xvg','w')
         cfVsCl.write("@    title \"Correlation Function\" \n")
         cfVsCl.write("@    xaxis  label \"Contour Length (s)(nm)\" \n")
         cfVsCl.write("@    yaxis  label \" C = <u(0).u(s)> \"\n")
         cfVsCl.write("@TYPE xy""\n")

         contourLen_tot_avg = avgUnequalArray(contourLen_tot)
         corrFunc_tot_avg = avgUnequalArray(corrFunc_tot)

         for i in range(len(contourLen_tot_avg)):
            cfVsCl.write(" %8.3f    %8.3f\n" %(contourLen_tot_avg[i][0], corrFunc_tot_avg[i][0]))
         lpTotalAllFrame += fit_exponential_decay(contourLen_tot_avg[:,0], corrFunc_tot_avg[:,0]) 

       keys = list(results.keys())
       nKeys = len(keys)
       keyLength = [(len(key)+3)*'=' for key in keys]
       output1 = open('lp_polymer.dat', 'w')
       output1.write(('    {}' + '\t\t{:>s}(nm)'*nKeys +'\n').format("Frame", *keys))
       output1.write(('    {:=^5}' + '\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       for i in range(len(results[keys[0]])):
           values = []
           for j in range(nKeys):
               values.append(results[keys[j]][i])
           output1.write(("    {:^5d}" + "\t\t{:.3f}"*nKeys + "\n").format(i+1, *values))

     
       # get the standard deviation of data   
       mean_of_keys = []
       std_of_keys= []
       for i in range(nKeys):
           avgKeys = mean(np.array(results[keys[i]]))
           stdKeys = standard_dev(np.array(results[keys[i]]))
           mean_of_keys.append(avgKeys)
           std_of_keys.append(stdKeys)
 
       # get merged list with average and sd in sequence
       merged_values = [j for i in zip(mean_of_keys, std_of_keys) for j in i] 
       output1.write(('    {:=^5}' + '\t\t{:>s}'*nKeys +'\n').format('', *keyLength))
       output1.write(('    {}' + '\t {:.3f} +/- {:.3f}'*nKeys +'\n').format("Avgerages: ", *merged_values))

       # get merged keys with corresponding values
       merged_all = []
       for i in range(len(keys)):
           merged_all.append(keys[i])
           for j in range(i*2, (i*2)+2):
               merged_all.append(merged_values[j])

       for i in range(len(merged_all)):
           if merged_all[i]=='rg':
              merged_all[i] = "Radius of Gyration"
           elif merged_all[i] == "ete":
              merged_all[i] = "End-to-End distance"
           elif merged_all[i] == "Lp_Exp":
              merged_all[i] = "Lp Exp. fit"

       print (("\n \t {:20s} : {:.3f} +/- {:.3f}"*nKeys +'\n').format(*merged_all))
       print("""\t Total Persistence Length = {:.3f}""".format(lpTotalAllFrame))
       output1.write("\t Total Persistence Length from all frame = {:.3f}".format(lpTotalAllFrame))

    end_time = timeit.default_timer()
    totalTimeTaken = end_time-start_time
    stdout.write("\t Time Taken = %i Seconds.\n \n" %totalTimeTaken)

def args():
    parser = ArgumentParser(description=__doc__, 
			   formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-f","--traj_file", help="""pass "multi" to use multiple .gro files for analysis.
                                                                   or.xtc(trajectory) file for analysis.""",
                        type=str)

    parser.add_argument('-g','--gro_file',
                        help="""
                        Single .gro file for analysis.
                        """,
                        type=str)

    parser.add_argument('-b','--begin_time',
		        help="""
                        Starting time for analysis.
                        \n        
			""",
                        type=int)
   
    parser.add_argument('-e','--end_time',
                        help="""
                        End time for analysis.
                        """,
                        type=int)

    parser.add_argument('-n','--index_no',
                        help="""
                        Index number for polymer chain.
                        """,
                        type=int)

    parser.add_argument('-s', '--topol_tpr',
                        help="""
                        topol.tpr file.
                        """,
                        type=str)

    parser.add_argument('-lp', help = """ 
                        Calculate persistence length of backbone polymer chain.Need to provide the 
                        backbone selection string. Eg. resname ST1 ST2 ST3 and name C31 C21 C11) .""", type= str)

    parser.add_argument('-lpk', help = """ 
                        Calculates segment persistence length along backbone of polymer chain 
                        using the projection of end to end vector into segment bond vector.Need to provide the 
                        backbone selection string. Eg. resname ST1 ST2 ST3 and name C31 C21 C11) .""", type= str)

    parser.add_argument('-rg', help = """ 
                        Calculate Rg of selected atoms (eg. resname ST1 ST2 and name C31 C21 C11) .""", type= str)

    parser.add_argument('-ee', nargs=2, help = """ 
                        Calculate End to End distance of polymer chain.Arguments needed are the 
                         name of end atoms of polymer. eg, "resname ST1 and name C11"  "resname ST3 and name C31" """, type= str)

    parser.add_argument('-eeacf', nargs = 2, help = """
                         Calculates the end to end vector auto correlation. arguments needed are the 
                         name of end atoms of polymer. eg, "resname ST1 and name C11"  "resname ST3 and name C31" """, type = str)
    parser.add_argument('-segacf', nargs = 2, help = """
                         Calculates the segmental vector auto correlation. Averaging is carried out
                         for vectors along the chain backbone separated by "sep" and at different starting times.
                         Arguments needed are the atoms in the backbone and int(sep), separation between atoms.
                         eg, "resname ST1 ST2 ST3 and name C31 C21 C11" 5 """, type = str)
    parser.add_argument('-tres', required='-segacf' in argv,
                        help = """ time in picoseconds for the segmental autocorrelation function to restart""", type=int) 

    parser.add_argument('-dt', required='-segacf' in argv,
                        help = """ separation time of frames in picoseconds in trajectory(traj.xtc). You should know it beforehand""", type= int)                        

    return parser, parser.parse_args() 

if __name__ =='__main__':
   parser, arg = args()
   main()
