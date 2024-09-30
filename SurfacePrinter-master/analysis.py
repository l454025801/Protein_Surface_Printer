#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:32:00 2023

@author: leon
"""
import numpy as np
from mathtools import *
import csv 

def count_and_print_domains(is_surface, protein, center, box, xyz, fileQp, fileQn, fileP, fileH):
    '''Count the number of protein surface domains: Positively charged, Negatively charged, Polar charge-neutral, Nonpolar'''
    '''The documents with spherical coordinates could be plotted using Wolfram Mathematica. '''
    '''The documents with cartesian coordinates could be visualized using commane (vmd protein.pdb surface.xyz).'''
    
    number_surface = np.zeros(5, dtype=int)      
    for ii in range(len(protein)):
            if is_surface[ii] == 0:     #non-surface atoms
                continue
            
            '''calcualte spherical coordinates'''
            if xyz==False: 
                theta, phi, dist = convert_to_spherical_coordinate(protein[ii].position, center, box)
            
            resname = protein[ii].resname
            if ( (resname =='ALA') or  #ALA    nonpolar amino acids, based on VMD
                 (resname =='LEU') or  #LEU
                 (resname =='VAL') or  #VAL
                 (resname =='ILE') or  #ILE
                 (resname =='PRO') or  #PRO
                 (resname =='PHE') or  #PHE
                 (resname =='MET') or  #MET
                 (resname =='TRP') or
                 (resname =='GLY')):   #TRP
                number_surface[3] += 1
                if xyz:
                    fileH.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileH.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            elif ((resname =='ARG') or     #ARG   positively charged
                 (resname =='HSP') or     #charged HIS = HSP
                 (resname =='HIS') or
                 (resname =='LYS')):     #LYS
                number_surface[0] += 1
                if xyz:
                    fileQp.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                 protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileQp.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            elif ((resname =='ASP') or  #ASP   negatively charged
                 (resname =='GLU')):  #GLU
                number_surface[1] += 1
                if xyz:
                    fileQn.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                 protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileQn.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            elif ((resname =='SER') or  #SER   polar neutral
                 (resname =='THR') or  #THR
                 (resname =='TYR') or  #VAL
                 (resname =='ASN') or  #ILE
                 (resname =='GLN') or
                 (resname =='HID') or 
                 (resname == 'HIE') or
                 (resname == 'CYS')):  #PRO:
                number_surface[2] += 1
                if xyz:
                    fileP.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileP.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            else:
                print(resname)
                number_surface[4] += 1
                
    return number_surface, fileQp, fileQn, fileP, fileH

def print_domain_probability(time, number_domains, file_format="csv"):
    '''print the probability of surface domains'''
    
    '''calcualte the probability (in percent) for each snapshot (row)'''
    if len(number_domains) == 1:
        total = np.sum(number_domains)/100.0
    else:
        total = np.sum(number_domains, axis=1)/100.0
        total = total.reshape(-1, 1)
    number_domains_prob = number_domains/total
    '''Average over all the snapshots (column)'''
    prob = np.mean(number_domains_prob, axis=0)
    std  = np.std(number_domains_prob, axis=0)
    print("***************************************************")
    print("The probability (in percent) is:")
    print("    positively charged: {:6.3f} +- {:.3f}".format(prob[0], std[0]))
    print("    negatively charged: {:6.3f} +- {:.3f}".format(prob[1], std[1]))
    print("    polar neutral     : {:6.3f} +- {:.3f}".format(prob[2], std[2]))
    print("    nonpolar          : {:6.3f} +- {:.3f}".format(prob[3], std[3]))
    print("    Unkown residues   : {:6.3f} +- {:.3f}".format(prob[4], std[4]))
    print("***************************************************")
    print("The number of protein surface atoms as a function of time (ps):")
    print("Time (ps)   Qp     Qn     P      H      Unkown" )
    if time == 0:
        for ii in range(len(number_domains)):
            print ("{:6g} {:6g} {:6g} {:6g} {:6g} {:6g}".format(time, *number_domains[ii], sep = ' '))
    else:
        for ii in range(len(number_domains)):
            if int(time[ii]%1000) == 0:    #print every 1ns
                print("{:6g} {:6g} {:6g} {:6g} {:6g} {:6g}".format(time[ii], *number_domains[ii], sep = ' '))
    print("***************************************************")
        
    '''print the output file'''
    if file_format == "dat":
        output_file = open("domain_distribution.dat","w")
        output_file.write("#Probability (in percent):\n")
        output_file.write("# Qp: {:6.3f}+-{:.3f}, Qn: {:6.3f}+-{:.3f}, P: {:6.3f}+-{:.3f}, H: {:6.3f}+-{:.3f}, Unkown: {:6.3f}+-{:.3f}\n"\
                          .format(prob[0], std[0], prob[1], std[1], prob[2], std[2], prob[3], std[3], prob[4], std[4]))
        output_file.write("#Tim (ps)   Qp     Qn     P      H      Unkonw\n" )
        if time == 0:
            for ii in range(len(number_domains)):
                output_file.write(" {:6g} {:6g} {:6g} {:6g} {:6g} {:6g}\n".format(time, *number_domains[ii], sep = ' '))
        else:
            for ii in range(len(number_domains)):
                output_file.write(" {:6g} {:6g} {:6g} {:6g} {:6g} {:6g}\n".format(time[ii], *number_domains[ii], sep = ' '))
        output_file.close()
        
    elif file_format == "csv":
        fields = ['Time(ps)','Qp','Qn','P','H','unknown']
        all_data = np.concatenate([np.asarray(time).reshape((len(time),1)), number_domains],1)
        with open("domain_distribution.csv","w",newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(fields)
            csvwriter.writerows(all_data)
        csvfile.close()
    return
