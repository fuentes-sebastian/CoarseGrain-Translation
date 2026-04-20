# CG_graph
A script to convert atoms in a .mol2, .itp or .top protein file into Coarse Grained (CG) beads using graph theory.

User define the input file. The output is a .txt file `final_dictionary.txt`, which is a list of dictionaries; each dictionary will correspond to a CG bead and will link the all-atom id with their corrsponding CG bead. Warnings will be written in the `warnings.txt` file. An example of usage is given for the `1aki.mol2` in the [tests/](https://github.com/fuentes-sebastian/CoarseGrain-Translation/tree/main/tests) folder.


> [!IMPORTANT]
> The input file must have all the hydrogens (particularly to correctly identify the serine residue) and all the bonds.


## How it works

The script reads a reference dictionary with the graphs of the residues `res_dictionary.txt` and compares them wth the residues of the input file, to assign each atom to its corresponding reference graph node. Then, the script correlates this information with the coarse grain dictionary `CG_dictionary.txt` to assign the original atom ids to beads in the CG model.

Here is an example of a CG bead dictionary generated in the `final_dictionary.txt`:

`{'CG_id': 14, 'atom_ids': [91, 92, 93, 94], 'CG_type': 'P2', 'res_name': 'CYS', 'res_id': 6}` 

This indicates that the 14th CG bead will: contain to the atoms `[91,92,93,94]` of the input file, be a "P2" bead type of the residue "CYS", and had an id of 6 in the input file.

## Using modified residues

Modified residues should be added to both the `res_dictionary.txt` and the `CG_dictionary.txt`, both on the [dictionaries/](https://github.com/fuentes-sebastian/CoarseGrain-Translation/tree/main/dictionaries) folder. The script `mol2_to_dict.py` can be used to create a reference dictionary from a .mol2 file of the modified residue; this script will generate a `res_dict.txt` with the new dictionary and that new dictionary can be pasted in the `res_dictionary.txt` file.

On the otherhand, the coarse grain information of the modified residue must be manually written in the `CG_dictionary.txt` file. 

Here is an example of a CG residue in the `CG_dictionary.txt`: 

`{'res_name': 'ARG', 'nodes': [[1, 2, 3, 4], [5, 6, 7], [8, 9, 10, 11]], 'names': ['P2', 'SC3', 'SQ3p']}`

This indicates that ARG residues will contain 3 CG beads with names `['P2', 'SC3', 'SQ3p']` and those beads will contain the atoms `[[1, 2, 3, 4], [5, 6, 7], [8, 9, 10, 11]]` of the reference dicitonary in `res_dictionary.txt`. The coarse grain dictionary can be assigned with only some atoms in the reference; for instance, in this example hydrogens are ignored.


> [!IMPORTANT]
> The nodes in the coarse grain dictionary must coincide with the nodes in the reference as this is the link between the all-atom model and the CG model. And the residue name must coincide with its name in the input_file
>

The script outputs a warning file with all the warnings found and a dictionary list with information for each coarse grained bead identified.
