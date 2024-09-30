#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:25:44 2023

@author: leon

Simple math tools
"""
import os 
from math import sqrt, acos
import numpy as np

'''Difference of vectors: vec1 - vec2 '''
def diff(vec1, vec2):
    return [coord1 - coord2 for coord1, coord2 in zip(vec1, vec2)]

'''To reset vector based on the PBC box'''
def pbc_dx(vec1, vec2, box, DIM=3):
    dv = diff(vec1, vec2)    
    for i in range(DIM):
        if dv[i] < -box[i]/2:
            dv[i] += box[i]
        elif dv[i] > box[i]/2:
            dv[i] -= box[i]
    return dv

'''To define a series of color representing polarity and distance''' 
def ChooseColor(a,x): 
    if a == 'r':
        return (1,0,0,x)
    if a == 'b':
        return (0,0,1,x)
    if a == 'k':
        return (0,0,0,x)
    if a == 'm':
        return (0.5,0,0.5,x)

def convert_dat2XYZ():
    '''Convert .data to .xyz file: add two header lines'''
    
    infile  = open("domain_location_Qp.dat","r")
    outfile = open("domain_location_Qp.xyz","w")
    nrLines = sum(1 for line in infile)
    outfile.write("{:}\n".format(nrLines))
    outfile.write("{:}\n".format('Positively charged domains'))
    infile.seek(0)
    for line in infile:
        outfile.write(line)
    infile.close()
    outfile.close()
    os.remove("domain_location_Qp.dat")

    infile  = open("domain_location_Qn.dat","r")
    outfile = open("domain_location_Qn.xyz","w")
    nrLines = sum(1 for line in infile)
    outfile.write("{}\n".format(nrLines))
    outfile.write("{}\n".format('negatively charged domains'))
    infile.seek(0)
    for line in infile:
        outfile.write(line)
    infile.close()
    outfile.close()
    os.remove("domain_location_Qn.dat")

    infile  = open("domain_location_P.dat","r")
    outfile = open("domain_location_P.xyz","w")
    nrLines = sum(1 for line in infile)
    outfile.write("{}\n".format(nrLines))
    outfile.write("{}\n".format('Polar charge-neutral domains'))
    infile.seek(0)
    for line in infile:
        outfile.write(line)
    infile.close()
    outfile.close()
    os.remove("domain_location_P.dat") 

    infile  = open("domain_location_H.dat","r")
    outfile = open("domain_location_H.xyz","w")
    nrLines = sum(1 for line in infile)
    outfile.write("{}\n".format(nrLines))
    outfile.write("{}\n".format('Nonpolar domains'))
    infile.seek(0)
    for line in infile:
        outfile.write(line)
    infile.close()
    outfile.close()
    os.remove("domain_location_H.dat")
    
    return 

def convert_to_spherical_coordinate(coordinate, center, box):
    '''convert cartesian coordinate to spherical coordinate''' 
    '''See Fig. 1 in [B.Qiao et al., PNAS 2019, 116, 19274-19281.]'''
    RAD2DEG = 57.29578049 
    dist_vec = pbc_dx(coordinate, center, box)
    dist =  np.linalg.norm(dist_vec)
    dist_uv = dist_vec/dist
    cosZ = dist_uv[2]
    cosX = dist_uv[0]/sqrt(dist_uv[0]*dist_uv[0] + dist_uv[1]*dist_uv[1])
    theta = RAD2DEG * acos(cosZ)
    phi   = RAD2DEG * acos(cosX) 
    if dist_uv[1] < 0:
        phi = 360 - phi
            
    return  theta, phi, dist