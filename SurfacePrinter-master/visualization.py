#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:33:44 2023

@author: leon
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from mathtools import *

def plot_domain_probability(time, number_domains):
    '''plot the probability of surface domains'''
    sns.set_style("darkgrid")

    plt.xlabel("Time (ps)")
    plt.ylabel("Number")
    plt.plot(time, number_domains[:, 0], '-b', label='Positively charged')
    plt.plot(time, number_domains[:, 1], '-r', label='Negatively charged')
    plt.plot(time, number_domains[:, 2], '-m', label='Polar charge-neutral')
    plt.plot(time, number_domains[:, 3], '-k', label='Nonpolar')
    plt.legend()
    plt.show()
    #plt.savefig('domain_probability.pdf')
    return

def Projection(is_surface, protein, center, box, xyz):
    '''Create 2D projection using Mercator Projection if spherical coordinate is used. Create 3D plot if Cartesian is used.'''
    
    Surface_info = [] # List of all residues with their name, polarity and coordinations
    
    ''' Seperate list for residues respect to their polarities'''
    Qp_info = []
    Qn_info = []
    P_info = []
    H_info = []
    distance = []
    
    if xyz == False:    
        for ii in range(len(protein)):
            if is_surface[ii] == 0:    #non-surface atoms
                continue
            
            #Create a dictionary containing resname, polarity and coordinations for each surface atom
            info = {} 
            

            #calcualte spherical coordinates
            theta, phi, dist = convert_to_spherical_coordinate(protein[ii].position, center, box)
            distance.append(dist)
    
            resname = protein[ii].resname
            if ( (resname[0] =='A' and resname[1] =='L' and resname[2] =='A') or  #ALA    nonpolar, based on VMD
                 (resname[0] =='L' and resname[1] =='E' and resname[2] =='U') or  #LEU
                 (resname[0] =='V' and resname[1] =='A' and resname[2] =='L') or  #VAL
                 (resname[0] =='I' and resname[1] =='L' and resname[2] =='E') or  #ILE
                 (resname[0] =='P' and resname[1] =='R' and resname[2] =='O') or  #PRO
                 (resname[0] =='P' and resname[1] =='H' and resname[2] =='E') or  #PHE
                 (resname[0] =='M' and resname[1] =='E' and resname[2] =='T') or  #MET
                 (resname[0] =='T' and resname[1] =='R' and resname[2] =='P')):   #TRP
                info['resname'] = protein[ii].resname
                info['polarity'] = 3
                info['coordinate'] = [theta,phi,dist]
                H_info.append(info)
            elif((resname[0] =='A' and resname[1] =='R' and resname[2] =='G') or  #ARG   positively charged
                 (resname[0] =='H' and resname[1] =='S' and resname[2] =='P') or  #charged HIS = HSP
                 (resname[0] =='L' and resname[1] =='Y' and resname[2] =='S') ):  #LYS
                info['resname'] = protein[ii].resname
                info['polarity'] = 0
                info['coordinate'] = [theta,phi,dist]
                Qp_info.append(info)
            elif((resname[0] =='A' and resname[1] =='S' and resname[2] =='P') or  #ASP   negatively charged
                 (resname[0] =='G' and resname[1] =='L' and resname[2] =='U') ):  #GLU
                info['resname'] = protein[ii].resname
                info['polarity'] = 1
                info['coordinate'] = [theta,phi,dist]
                Qn_info.append(info)
            else:                                                  #polar charge-neutral
                info['resname'] = protein[ii].resname
                info['polarity'] = 2
                info['coordinate'] = [theta,phi,dist]
                P_info.append(info)
            Surface_info.append(info)
        
            '''print Surface_info'''
        dist_max = max(distance)   
        for residue in Surface_info:
            if residue['polarity'] == 0 :
                plt.plot(residue['coordinate'][0],residue['coordinate'][1],'o',color = ChooseColor('b',residue['coordinate'][2]/dist_max)) 
            if residue['polarity'] == 1 :
                plt.plot(residue['coordinate'][0],residue['coordinate'][1],'o',color = ChooseColor('r',residue['coordinate'][2]/dist_max))
            if residue['polarity'] == 2 :
                plt.plot(residue['coordinate'][0],residue['coordinate'][1],'o',color = ChooseColor('m',residue['coordinate'][2]/dist_max))
            if residue['polarity'] == 3 :
                plt.plot(residue['coordinate'][0],residue['coordinate'][1],'o',color = ChooseColor('k',residue['coordinate'][2]/dist_max))
        plt.xlabel("\u03B8 (deg.)")       #theta
        plt.ylabel("\u03C6 (deg.)")       #phi
        plt.xticks(np.arange(0, 181, 30))
        plt.yticks(np.arange(0, 361, 60))
    
        Qp_legend = mpatches.Circle((175,275),color = 'b', label ='Positively charged')
        Qn_legend = mpatches.Circle((175,275),color = 'r', label ='Negatively charged')
        P_legend  = mpatches.Circle((175,275),color = 'm', label ='Polar charge-neutral')
        H_legend  = mpatches.Circle((175,275),color = 'k', label ='Nonpolar')
        plt.legend(handles = [Qp_legend,Qn_legend,H_legend,P_legend],loc = 'upper left', bbox_to_anchor=(1.05, 0.5))
        
    else:
        for ii in range(len(protein)):
            if is_surface[ii] == 0:    #non-surface atoms
                continue
            
            #Create a dictionary containing resname, polarity and coordinations for each surface atom
            info = {} 
            x = protein[ii].position[0]
            y = protein[ii].position[1]
            z = protein[ii].position[2]
            
            resname = protein[ii].resname
            if ( (resname[0] =='A' and resname[1] =='L' and resname[2] =='A') or  #ALA    nonpolar, based on VMD
                 (resname[0] =='L' and resname[1] =='E' and resname[2] =='U') or  #LEU
                 (resname[0] =='V' and resname[1] =='A' and resname[2] =='L') or  #VAL
                 (resname[0] =='I' and resname[1] =='L' and resname[2] =='E') or  #ILE
                 (resname[0] =='P' and resname[1] =='R' and resname[2] =='O') or  #PRO
                 (resname[0] =='P' and resname[1] =='H' and resname[2] =='E') or  #PHE
                 (resname[0] =='M' and resname[1] =='E' and resname[2] =='T') or  #MET
                 (resname[0] =='T' and resname[1] =='R' and resname[2] =='P')):   #TRP
                info['resname'] = protein[ii].resname
                info['polarity'] = 3
                info['coordinate'] = [x,y,z]
                H_info.append(info)
            elif((resname[0] =='A' and resname[1] =='R' and resname[2] =='G') or  #ARG   positively charged
                 (resname[0] =='H' and resname[1] =='S' and resname[2] =='P') or  #charged HIS = HSP
                 (resname[0] =='L' and resname[1] =='Y' and resname[2] =='S') ):  #LYS
                info['resname'] = protein[ii].resname
                info['polarity'] = 0
                info['coordinate'] = [x,y,z]
                Qp_info.append(info)
            elif((resname[0] =='A' and resname[1] =='S' and resname[2] =='P') or  #ASP   negatively charged
                 (resname[0] =='G' and resname[1] =='L' and resname[2] =='U') ):  #GLU
                info['resname'] = protein[ii].resname
                info['polarity'] = 1
                info['coordinate'] = [x,y,z]
                Qn_info.append(info)
            else:                                                  #polar charge-neutral
                info['resname'] = protein[ii].resname
                info['polarity'] = 2
                info['coordinate'] = [x,y,z]
                P_info.append(info)
            Surface_info.append(info)
    
        '''print Surface_info'''
        ax = plt.axes(projection = '3d')
        for residue in Surface_info:
            if residue['polarity'] == 0 :
                ax.scatter3D(residue['coordinate'][0],residue['coordinate'][1],residue['coordinate'][2],color = 'b') 
            if residue['polarity'] == 1 :
                ax.scatter3D(residue['coordinate'][0],residue['coordinate'][1],residue['coordinate'][2],color = 'r')
            if residue['polarity'] == 2 :
                ax.scatter3D(residue['coordinate'][0],residue['coordinate'][1],residue['coordinate'][2],color = 'm')
            if residue['polarity'] == 3 :
                ax.scatter3D(residue['coordinate'][0],residue['coordinate'][1],residue['coordinate'][2],color = 'k')
        ax.set_xlabel('X (\u212B)', size = '15', labelpad = 10)
        ax.set_ylabel('Y (\u212B)', size = '15', labelpad = 10)
        ax.set_zlabel('Z (\u212B)', size = '15', labelpad = 10)
        ax.tick_params(labelsize = '15')
        ax.dist = 10
        
        Qp_legend = mpatches.Circle((175,275),color = 'b', label ='Positively charged')
        Qn_legend = mpatches.Circle((175,275),color = 'r', label ='Negatively charged')
        P_legend  = mpatches.Circle((175,275),color = 'm', label ='Polar charge-neutral')
        H_legend  = mpatches.Circle((175,275),color = 'k', label ='Nonpolar')
        ax.legend(handles = [Qp_legend,Qn_legend,H_legend,P_legend],loc = 'upper left', bbox_to_anchor=(1.05, 0.5))
        return
    