#!/usr/bin/env python
__description__ = '''To print the protein surface domains '''
__author__ = "Baofu Qiao, Yang Li"
__date__  = "031219"
__usage__ = "python3 protein_surface_domain.py -trj trajectory.xtc -top topology.tpr -[no]sphere"
__version__ = "1"

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
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (6, 5),
         'axes.labelsize': 30,
         'axes.titlesize':30,
         'xtick.labelsize':20,
         'ytick.labelsize':20}
pylab.rcParams.update(params)

RAD2DEG = 57.29578049 

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
                                    help = "Align protein based on index group 1 (default).")
        
    args = parser.parse_args()
    return args
    
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
    

def count_and_print_domains(is_surface, protein, center, box, xyz, fileQp, fileQn, fileP, fileH):
    '''Count the number of protein surface domains: Positively charged, Negatively charged, Polar charge-neutral, Nonpolar'''
    '''The documents with spherical coordinates could be plotted using Wolfram Mathematica. '''
    '''The documents with cartesian coordinates could be visualized using commane (vmd protein.pdb surface.xyz).'''
    
    number_surface = np.zeros(4, dtype=int)      
    for ii in range(len(protein)):
            if is_surface[ii] == 0:     #non-surface atoms
                continue
            
            '''calcualte spherical coordinates'''
            if xyz==False: 
                theta, phi, dist = convert_to_spherical_coordinate(protein[ii].position, center, box)
            
            resname = protein[ii].resname
            if ( (resname[0] =='A' and resname[1] =='L' and resname[2] =='A') or  #ALA    nonpolar amino acids, based on VMD
                 (resname[0] =='L' and resname[1] =='E' and resname[2] =='U') or  #LEU
                 (resname[0] =='V' and resname[1] =='A' and resname[2] =='L') or  #VAL
                 (resname[0] =='I' and resname[1] =='L' and resname[2] =='E') or  #ILE
                 (resname[0] =='P' and resname[1] =='R' and resname[2] =='O') or  #PRO
                 (resname[0] =='P' and resname[1] =='H' and resname[2] =='E') or  #PHE
                 (resname[0] =='M' and resname[1] =='E' and resname[2] =='T') or  #MET
                 (resname[0] =='T' and resname[1] =='R' and resname[2] =='P')):   #TRP
                number_surface[3] += 1
                if xyz:
                    fileH.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileH.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            elif((resname[0] =='A' and resname[1] =='R' and resname[2] =='G') or     #ARG   positively charged
                 (resname[0] =='H' and resname[1] =='S' and resname[2] =='P') or     #charged HIS = HSP
                 (resname[0] =='L' and resname[1] =='Y' and resname[2] =='S') ):     #LYS
                number_surface[0] += 1
                if xyz:
                    fileQp.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                 protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileQp.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            elif((resname[0] =='A' and resname[1] =='S' and resname[2] =='P') or  #ASP   negatively charged
                 (resname[0] =='G' and resname[1] =='L' and resname[2] =='U') ):  #GLU
                number_surface[1] += 1
                if xyz:
                    fileQn.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                 protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileQn.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            else:                                                  #polar neutral
                number_surface[2] += 1
                if xyz:
                    fileP.write("{:4s}  {:.3f}  {:.3f}  {:.3f}\n".format(protein[ii].name, \
                                protein[ii].position[0], protein[ii].position[1], protein[ii].position[2]))
                else:
                    fileP.write("{:6g}  {:6g}\n".format(int(theta), int(phi)))
            
    return number_surface, fileQp, fileQn, fileP, fileH
    

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
        


def convert_to_spherical_coordinate(coordinate, center, box):
    '''convert cartesian coordinate to spherical coordinate''' 
    '''See Fig. 1 in [B.Qiao et al., PNAS 2019, 116, 19274-19281.]'''
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


def print_domain_probability(time, number_domains):
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
    print("***************************************************")
    print("The number of protein surface atoms as a function of time (ps):")
    print("Time (ps)   Qp     Qn     P      H" )
    if time == 0:
        for ii in range(len(number_domains)):
            print ("{:6g} {:6g} {:6g} {:6g} {:6g}".format(time, *number_domains[ii], sep = ' '))
    else:
        for ii in range(len(number_domains)):
            if int(time[ii]%1000) == 0:    #print every 1ns
                print("{:6g} {:6g} {:6g} {:6g} {:6g}".format(time[ii], *number_domains[ii], sep = ' '))
    print("***************************************************")
        
    '''print the output file'''
    output_file = open("domain_distribution.dat","w")
    output_file.write("#Probability (in percent):\n")
    output_file.write("# Qp: {:6.3f}+-{:.3f}, Qn: {:6.3f}+-{:.3f}, P: {:6.3f}+-{:.3f}, H: {:6.3f}+-{:.3f}\n".format(
                               prob[0], std[0], prob[1], std[1], prob[2], std[2], prob[3], std[3]))
    output_file.write("#Tim (ps)   Qp     Qn     P      H\n" )
    if time == 0:
        for ii in range(len(number_domains)):
            output_file.write(" {:6g} {:6g} {:6g} {:6g} {:6g}\n".format(time, *number_domains[ii], sep = ' '))
    else:
        for ii in range(len(number_domains)):
            output_file.write(" {:6g} {:6g} {:6g} {:6g} {:6g}\n".format(time[ii], *number_domains[ii], sep = ' '))
    output_file.close()

    
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

    
if __name__ == "__main__":
    #parse command line
    args = ParseOptions()
    #do the calculations
    do_calculation(args)