# CG_graph
A script to convert atoms in a .mol2, .itp or .top protein file into Coarse Grained (CG) beads using graph theory. 

The input must have all hydrogens (particularly to correctly identify the serine residue) and bonds of the protein atoms. 
For the moment it does not accept modified residues

The script takes an input file (.mol2, .itp or .top), uses a reference dictionary of the residues and a CG dictionary of the CG beads for each residue. 

The output is a warning file with all the warnings found and a dictionary list for each coarse grained bead.
