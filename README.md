# CG_graph
A script to convert atoms in a .mol2, .itp or .top protein file into Coarse Grained (CG) beads using graph theory. 

The input file must have all the hydrogens (particularly to correctly identify the serine residue) and bonds of the protein atoms. User define the input and the output is a .txt file "final_dictionary.txt" which is a list of dictionaries of each bead.

The script reads a reference dictionary with the graphs of the residues (res_dictionary.txt) and compares each residue of the input file with the reference graphs, to assign each atom to its corresponding reference graph node. Then, the script correlates this information with the CG dictionary (CG_dictionary.txt) to assign the original atom ids to beads in the CG model

Modified residues can be added to both the res_dictionary and the CG_dictionary.

The script outputs a warning file with all the warnings found and a dictionary list with information for each coarse grained bead identified.
