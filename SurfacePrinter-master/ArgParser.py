#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:28:11 2023

@author: leon

Functions used to parse arguments from user input
"""
import argparse

def ParseOptions():
    '''Parse the command and return a set of options'''
    parser = argparse.ArgumentParser(description="Print protein surface domains. " 
                                     "(1) the numbers of different kinds of surface domains in domain_distribution.dat. " 
                                     "(2) Location (in XYZ of spherical coordinate) of all kinds of domains: domain_location.xyz (.dat). "
                                     "Use [vmd protein.pdb domain_location.xyz] to visualize the 3D distribution." ,
                                                 usage='use "%(prog)s --help" for more information')
    parser.add_argument("-trj",     dest="trj_filename", default='./Gromacs/aligned_trj.xtc',  
                                    action="store",  help="read a trajectory file ")
    parser.add_argument("-top",     dest="top_filename", default='./Gromacs/T298-3.gro',  
                                    action="store",  help="read a topology file")
    parser.add_argument("-b",       dest="start_time", default=0, type=int, 
                                    action="store",  help="First time to Read trajectory")
    parser.add_argument("-e",       dest="end_time", default=float('2000'), type=int,  #default=float('inf')
                                    action="store",  help="Last time to Read trajectory")
    parser.add_argument("-skip",    dest="skip", default=10, type=int, 
                                    action="store", help="how often to read frame (default: 1)")
    
    parser.add_argument("-ndx1",    dest="index_group1", default='backbone and name N CA C', type=str,   
                                    action="store",  help="index group 1")
    parser.add_argument("-ndx2",    dest="index_group2", default='name OW NA CL and around 5 protein', type=str,  
                                    action="store",  help="index group 2")
    
    parser.add_argument("-xyz",     dest="xyz", action="store_true", default=True, 
                                    help="print domain distributio in cartesian coordinate (default), or spherical coordinate")
    parser.add_argument("-com",     dest="COM", action="store_true", default=True,
                                    help="use center-of-mass (default), or geometric center")
    parser.add_argument("-align",   dest="alignment", action = "store_true", default=False, 
                                    help = "(True/False) Align protein based on index group 1. Default=False")
        
    args = parser.parse_args()
    return args
