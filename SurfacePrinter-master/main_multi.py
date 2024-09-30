#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:30:57 2023

@author: leon
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from math import sqrt, acos

'''----- Read MDAnalysis -------------'''
import MDAnalysis
from MDAnalysis import transformations
from MDAnalysis.analysis import align

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pylab as pylab

from mathtools import *
from ArgParser import *
from analysis import *
from visualization import *

def do_calculation(args):
    '''Read the input parameters'''
    trj_filename = args.trj_filename
    top_filename = args.top_filename

    index_group1 = args.index_group1
    index_group2 = args.index_group2      

    start_time = args.start_time
    end_time   = args.end_time
    skip = args.skip
    COM = args.COM
    xyz = args.xyz
    
    alignment = args.alignment
    
    print ('trajectory file: %s' % trj_filename)
    print ('  Topology file: %s' % top_filename)
    
    if top_filename == 'null':
        # for a single frame
        print("*******************************************************************************")
        print("The following parameters are used:") 
        print("  Single Frame")
        print("  protein center-of-mass: {}; Cartesian coordinate: {};".format(COM, xyz))
        print("  Index group 1: {}".format(index_group1))
        print("        group 2: {}".format(index_group2))
        print("*******************************************************************************")
        
        u = MDAnalysis.Universe(trj_filename)
        protein = u.select_atoms(index_group1)
        neighbors = u.select_atoms(index_group2)
        box = u.dimensions
        
        if (COM == 1):
            center = protein.center_of_mass()
        else:
            center = protein.center_of_geometry()
            
        #check if atom correctly select
        group1 = u.select_atoms(index_group1)
        group2 = u.select_atoms(index_group2)
        
        if len(group1) == 0:
            print ('Empty atom selection: Group1; Check top file for atom name')
            sys.exit()
        if len(group2) == 0:
            print ('Empty atom selection: Group2; Check top file for atom name')
            sys.exit()
        
        print ('-------------------------------Start Calculation-------------------------------')
        
        fileQp = open("domain_location_Qp.dat","w")
        fileQn = open("domain_location_Qn.dat","w")
        fileP  = open("domain_location_P.dat", "w")
        fileH  = open("domain_location_H.dat", "w")
        
        clust_index = []
        clust_size = []
        for i in range(len(neighbors)):
            clust_index.append(i)
            clust_size.append(1)
        # loop over every atoms, if they are adjacent, put them in a cluster
        print ("number of selected group2 atoms: %i" % len(neighbors))
        for i in range(len(neighbors)):
            cluster_i = clust_index[i]
            for j in range(i+1, len(neighbors)):
                cluster_j = clust_index[j]
                if cluster_i != cluster_j:
                    dist_vec = pbc_dx(neighbors[j].position, neighbors[i].position, box)
                    dist = np.linalg.norm(dist_vec)
                    if dist <= 3.5:   #Donor-acceptor distance for Hydrogen Bonds between water
                        for k in range(len(neighbors)):
                            if clust_index[k] == cluster_j:
                                clust_size[cluster_j] -= 1
                                clust_index[k] = cluster_i
                                clust_size[cluster_i] += 1
        # Create a list for neighbors excluding atoms inside the molecule
        new_neighbors = [] 
        for i in range(len(neighbors)):
            cluster = max(set(clust_index), key = clust_index.count)
            if clust_index[i] == cluster:
                new_neighbors.append(neighbors[i])

        print ("Finish neighbor atoms calculations")
        is_surface = np.zeros(len(protein), dtype=int)
        for ii in range(len(new_neighbors)):
            dist_vec = protein.positions - new_neighbors[ii].position.reshape(1, -1)
            dist = [np.linalg.norm(dv) for dv in dist_vec]
            #index of the protein backbone atom with smallest distance to the selected neighbor ii
            ndx = np.argmin(dist) 
            #label the protein backbone atom as a surface atom
            is_surface[ndx] = 1 

        Projection(is_surface, protein, center, box, xyz)
        
        '''count the number of protein surface atoms: num_curr;
         And print the .xyz files, which can be read by VMD after convert_dat2xyz()'''
        number_curr, fileQp, fileQn, fileP, fileH = count_and_print_domains(is_surface, protein, center, box, xyz, fileQp, fileQn, fileP, fileH)
        #update the numbers of different types of domains: number_domains
        number_domains = np.vstack([number_curr])
        
        fileQp.close()
        fileQn.close()
        fileP.close()
        fileH.close()
        
    
        '''convert .dat (fileQp/Qn/P/H) to .xyz file. Use the command below to visualize the distribution of the surface domains
        [vmd protein.pdb domain_location_Qp/Qn/P/H.xyz]'''
        if xyz:
            convert_dat2XYZ()
    
        print_domain_probability(0, number_domains)   
        
    else:
        # for multiframes
        print("*******************************************************************************")
        print("The following parameters are used:") 
        print("  Multiframes")
        print("  start_time={} ns, end_time={} ns, skip_frame={}".format(start_time/1000, end_time/1000, skip))
        print("  protein center-of-mass: {}; Cartesian coordinate: {}; Protein alignment: {}".format(COM, xyz, alignment))
        print("  Index group 1: {}".format(index_group1))
        print("        group 2: {}".format(index_group2))
        print("*******************************************************************************")
            
        '''Align the protein based on the protein backbone atoms'''
        if alignment:
            print ('Starting alignment')
            print ('Might be TIME-CONSUMING; Expect to take minutes depending on your trajectory length')
            
            mobile = MDAnalysis.Universe(top_filename, trj_filename)
            system = mobile.select_atoms('all')
            protein = mobile.select_atoms('protein')
            not_protein = mobile.select_atoms('not protein')        
                    
            #unwrap protein to make it a whole molecule, center protein, and wrap the whole system
            total_steps = len(mobile.trajectory)
            total_time = mobile.trajectory[-1].time
            time_step = (total_steps-1)/total_time
            
            print ('Number of total frames (initial): %i' % total_steps)
            print ('            Total time (initial): %f (ps)' % total_time)
            print ('                       Time step: %f' % time_step)
            
            with MDAnalysis.Writer("centered_traj.xtc", system.n_atoms) as W:
                for ts in mobile.trajectory[int(start_time*time_step):int(end_time*time_step+2):1]:
                    protein.unwrap(compound='fragments')
                    ts = MDAnalysis.transformations.center_in_box(protein, center='mass')(ts)
                    not_protein.wrap(compound='fragments')
                    system = mobile.select_atoms('all')
                    W.write(system)
            
            #Save the 1st frame, which is already centered, as the reference for protein alignment
            mobile_2 = MDAnalysis.Universe(top_filename, "centered_traj.xtc")
            system = mobile_2.select_atoms('all')
            with MDAnalysis.Writer("reference.pdb") as pdb:
                pdb.write(system)
            ref = MDAnalysis.Universe(top_filename, 'reference.pdb')
            #Align protein
            aligned_trj_filename = 'aligned_trj.xtc'
            align.AlignTraj(mobile_2, ref, select = index_group1, filename = aligned_trj_filename).run()    
            u = MDAnalysis.Universe(top_filename, aligned_trj_filename)
            os.remove("centered_traj.xtc")
            print ("-----------------------------Successfully Aligned !----------------------------")
            
        else:
            u = MDAnalysis.Universe(top_filename, trj_filename)
        

        print ('  Number of total frames: %i' % len(u.trajectory))
        print ('Total time of trajectory: %f (ps)' % u.trajectory[-1].time)
        
        #check if atom correctly select
        group1 = u.select_atoms(index_group1)
        group2 = u.select_atoms(index_group2)
        
        if len(group1) == 0:
            print ('Empty atom selection: Group1; Check top file for atom name')
            sys.exit()
        if len(group2) == 0:
            print ('Empty atom selection: Group2; Check top file for atom name')
            sys.exit()
    
        print ('-------------------------------Start Calculation-------------------------------')
    
        '''Files to save the protein surface domain locations'''
        fileQp = open("domain_location_Qp.dat","w")
        fileQn = open("domain_location_Qn.dat","w")
        fileP  = open("domain_location_P.dat", "w")
        fileH  = open("domain_location_H.dat", "w")
        
        time, frame = [], []
        for ts in u.trajectory:
            '''get time and frame, and control which frames to read'''
            curr_time = u.trajectory.time
            curr_frame = u.trajectory.frame
            if curr_time < start_time or curr_frame%skip != 0:
                continue
            elif curr_time > end_time:
                break
            time.append(curr_time)
            frame.append(curr_frame)
    
            box = u.dimensions
           
            '''define protein (backbone atoms) and calculate its center'''
            protein = u.select_atoms(index_group1)
            if (COM == 1):
                center = protein.center_of_mass()
            else:
                center = protein.center_of_geometry()
            
            '''Define the neighbors (mostly water)'''
            neighbors = u.select_atoms(index_group2)
            #refine the neighbors to exclude those inside protein (a small amount of water goes inside proteins)
            clust_index = []
            clust_size = []
            for i in range(len(neighbors)):
                clust_index.append(i)
                clust_size.append(1)
            # loop over every atoms, if they are adjacent, put them in a cluster
            for i in range(len(neighbors)):
                cluster_i = clust_index[i]
                for j in range(i+1, len(neighbors)):
                    cluster_j = clust_index[j]
                    if cluster_i != cluster_j:
                        dist_vec = pbc_dx(neighbors[j].position, neighbors[i].position, box)
                        dist = np.linalg.norm(dist_vec)
                        if dist <= 3.5:   #Donor-acceptor distance for Hydrogen Bonds between water
                            for k in range(len(neighbors)):
                                if clust_index[k] == cluster_j:
                                    clust_size[cluster_j] -= 1
                                    clust_index[k] = cluster_i
                                    clust_size[cluster_i] += 1
            # Create a list for neighbors excluding atoms inside the molecule
            new_neighbors = [] 
            for i in range(len(neighbors)):
                cluster = max(set(clust_index), key = clust_index.count)
                if clust_index[i] == cluster:
                    new_neighbors.append(neighbors[i])
             
            '''Print some info'''           
            if curr_frame%100 == 0:
                print("Protein backbone center: ({:.3f}, {:.3f}, {:.3f}) Angstrom at time {} ps".format(center[0], center[1], center[2], curr_time))
                print("    Found {} neighbors in group ({})".format(len(new_neighbors), index_group2))
    
            '''define the protein surface atoms'''
            is_surface = np.zeros(len(protein), dtype=int)
            for ii in range(len(new_neighbors)):
                dist_vec = protein.positions - new_neighbors[ii].position.reshape(1, -1)
                dist = [np.linalg.norm(dv) for dv in dist_vec]
                #index of the protein backbone atom with smallest distance to the selected neighbor ii
                ndx = np.argmin(dist) 
                #label the protein backbone atom as a surface atom
                is_surface[ndx] = 1 
    
            '''Plot coordinate'''
            Projection(is_surface, protein, center, box, xyz)
            
            '''count the number of protein surface atoms: num_curr;
             And print the .xyz files, which can be read by VMD after convert_dat2xyz()'''
            number_curr, fileQp, fileQn, fileP, fileH = count_and_print_domains(is_surface, protein, center, box, xyz, fileQp, fileQn, fileP, fileH)
            #update the numbers of different types of domains: number_domains
            if curr_frame == 0:
                number_domains = number_curr
            else:
                number_domains = np.vstack([number_domains, number_curr])
        '''Done reading the trajectory'''
        print('{:d} frames found! \n'.format(len(frame)))
    

        #plt.savefig('domain_distribution.pdf')    
        fileQp.close()
        fileQn.close()
        fileP.close()
        fileH.close()
        
        '''convert .dat (fileQp/Qn/P/H) to .xyz file. Use the command below to visualize the distribution of the surface domains
        [vmd protein.pdb domain_location_Qp/Qn/P/H.xyz]'''
        if xyz:
            convert_dat2XYZ()
        
        '''print and plot distribution probability of surface domains '''
        plt.show()
        print_domain_probability(time, number_domains)   
        plot_domain_probability(time, number_domains)
        plt.show()
        
    return 

if __name__ == "__main__":
    #parse command line
    args = ParseOptions()
    #do the calculations
    do_calculation(args)