#! /usr/bin/env python

""" 
   Description:
   ===========

   Calculates the survival probability as a function of time and the rate of 
   folding in units of 1/ns. A protein is said to be in unfolded state if the 
   fraction of native contacts in a given configuration, Q < 0.6, otherwise folded.

   Usage: ./survival_probability -f 1jw2_traj -o <output_filename> -p  
"""

__autohor = "Hari Sharma <hari.sharma@uconn.edu>"

from math import sqrt
import numpy as np
import timeit
from os import path
import os
from sys import argv, stdout
import matplotlib.pyplot as plt
from pylab import polyfit
from scipy.optimize import curve_fit
import argparse
from argparse import ArgumentParser, RawDescriptionHelpFormatter, FileType

np.set_printoptions(threshold=np.nan)
start_time = timeit.default_timer()

def weibull(x, a):
    return np.exp(-x**a)

def fit_weibull(x, y):
    a = curve_fit(weibull, x, y)[0][0]
    return a
   

def expfunc(x, a):
    return np.exp(-a*x)

def fit_exponential_decay(x, y):
    r"""Fit a function to an exponential decay
    .. math::  y = \exp(-a*x)

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
    a = curve_fit(expfunc, x, y)[0][0]
    return a

def read_trajectory(filename):
    """
    Reads the trajectory and converts the Q values in probabilistic 
    interpretation as a function of time. 
    i.e., if Q < 0.6, probability = 1; remains in unfolded state. 
    """ 
    count = 0
    time  = [0.0]
    prob = []

    while os.path.exists(filename + "%s_Q.dat" %count):
       file = filename + str(count) + "_Q.dat"
       f = open(file, 'r')
       lines = f.readlines()
       f.close()
       nlines = sum(1 for line in lines)

       p = np.zeros(nlines+1)
       p[0] += 1
    
       if count == 0:
          for i in range(nlines):
              t = [int(float(x)) for x in lines[i].split()][0]
              time.append(t)

       for i in range(nlines):
              q = [float(x) for x in lines[i].split()][1]
              if q < 0.6:
                 p[i+1] += 1
              else:
                 break

       p = np.array(p)
       prob.append(p)
       count += 1

   # convert integration time steps into real time (15 fs = 0.015 ns)  
    time = np.array(time)*0.000015
    
    prob = np.array(prob)    
    prob_normalized = np.mean(prob, axis=0)

   # truncate and return only values of probability greater than 0 and corresponding
   #  time 
    time_truncated = []
    prob_truncated = []
    
    for i in range(len(time)):
        if prob_normalized[i] > 0.0:
           time_truncated.append(time[i])
           prob_truncated.append(prob_normalized[i])
        else:
           break      

    return(np.array(time_truncated), np.array(prob_truncated))

def write_survival_prob(time, prob):
    """
       Writes the result of survival probability in a .xvg file 
    """
    sp = open(arg.output, 'w')
    sp.write("@    title \"Survival Probability vs time\" \n")
    sp.write("@    xaxis  label \"t(ns) \" \n")
    sp.write("@    yaxis  label \" Survival Probability p(t) \"\n")
    sp.write("@TYPE xy""\n")
    for i in range(len(time)):
        sp.write(" %8.3f    %8.3f\n" %(time[i], prob[i]))
     

def main():
   if len(argv) <2:
      print(parser.print_help())
      exit(1)

   if arg.input:
      filename = arg.input 
      time, prob = read_trajectory(filename)
      decay_rate = fit_exponential_decay(time, prob)
      stdout.write("\n Decay rate of protein folding = %.3f ns^-1 \n"%decay_rate)
 
   if arg.output:
      write_survival_prob(time, prob) 
      stdout.write("\n Survival probability vs time(ns) written in " +  arg.output + " file. \n\n")
  
   if arg.plot:
      yy = expfunc(time, decay_rate)
      #weibull fit
      wb = fit_weibull(time, prob)  
      zz = weibull(time, wb)

      plt.plot(time, prob,'bo', label="probability")
      plt.plot(time, yy, label="exp. fit (exp(-x*t))", color="red", linewidth=2.5)
      plt.plot(time, zz, label="weibull fit (exp(-x^t))", color="black", linewidth=2.5)

      plt.xlabel('time (ns)')
      plt.ylabel('Survival Probability P(t)')
      plt.title('Survival Probability')
      plt.legend()
      plt.show()


def args():
    parser = ArgumentParser(description=__doc__,
                           formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument('-f','--input',
                        help="""
                        name of trajectory files as input
                        """,
                        type=str)

    parser.add_argument('-o','--output',
                        help="""
                        outputs the result in filename you provided eg. sp.xvg
                        """,type=str)

    parser.add_argument('-p','--plot', action='store_true',
                        help="""
                           plots the result of survival probability as a function of time.
                        """)
    return parser, parser.parse_args()


if __name__ =='__main__':
   parser, arg = args()
   main()

