Protein Surface Printer

A versatile program that analyzes protein surfaces from MD trajectory.

Usage:


Python /SurfacePrinter/main.py -trj TRAJ_FILE -top TOPO_FILE -ndx1 SELECTED_GROUP -ndx2 REF_GROUP -b BEGIN_TIME -e END_TIME -skip FREQUENCY


the -ndx1 & -ndx2 are strings that used to select atoms following MDAnalysis syntax
https://docs.mdanalysis.org/stable/documentation_pages/selections.html

The program will analyze the surface of molecules in -ndx1 that is in contact with atoms in -ndx2.

Ex:
-ndx1 "protein" -ndx2 "water" will analyze the solvent accessible surface area of the protein in water.

For more details, please refer to the article:
J. Chem. Inf. Model. 2020, 60, 10, 5255–5264, DOI: doi.org/10.1021/acs.jcim.0c00582
