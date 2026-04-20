#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 03:32:49 2026

@author: sebas_fu
"""

# takes a mol2 file and creates a dictionary that can be converted to a graph

# ___________________________ Variables ___________________________ 
input_file='molecule.mol2' # file with mol2 format
res_name='RES' # name of residue, defined by user

# parsing
weight_dict={'H':1, 'C':12, 'N':14, 'O':16, 'S':32} # needed because mol2 files dont have atoms weight

def parse_atoms_mol2(line):
    atom=line.split()
    element=atom[1][0] # uses first letter of the atom name

    # sometimes atom name starts with a number, 
    # in that case, it will use first letter of the atom type
    if element not in weight_dict.keys():
        element=atom[5][0]
    
    # if element is not in either column:
    if element not in weight_dict.keys():        
        warn=f'atom with id {int(atom[0])} has element {element} which is not in dictionay of weights, wheight = 0 applied'
        print(warn)
        warnings.append(warn)
        weight=0
    if element in weight_dict.keys():
        weight=weight_dict[element]

    atom_info = {
        'id': int(atom[0]),
        'atom_name':atom[1],
        'element': element,
        'res_id': int(atom[6]),
        'res_name': atom[7] if len(atom) > 6 else 'Unk',
        'weight':weight}
    return(atom_info)

def bonds_prase_mol2(line):
    bond=line.split()
    bond_info = {
        'id': int(bond[0]),
        'a_i': int(bond[1]),
        'a_j': int(bond[2]),
        'bond_type' : int(bond[3])}
    return(bond_info)

# ___________________________ Parsing ___________________________ 

with open(f"{input_file}", 'r') as f:
    lines = f.readlines()

warnings=[] # warning list
dict_list=[] # list to append dictionaries of residues
atom_data=[]
bond_data=[]
section=None
for line in lines:
    line = line.strip() # eliminate space between lines
    if line.startswith("@<TRIPOS>"): # this ignores section line
        section = line
        continue

    # Atoms data prasing
    if section == "@<TRIPOS>ATOM":
        a_info=parse_atoms_mol2(line)
        atom_data.append(a_info)

    if section == "@<TRIPOS>BOND":
        bond=line.split()
        b_info=bonds_prase_mol2(line)
        bond_data.append(b_info)

# lists to write in dictionary
bond_list=[]
for b in bond_data:
    bond_list.append((b['a_i'],b['a_j']))

element_list=[]
node_list=[]
weight_list=[]
atom_names=[]
for a in atom_data:
    element_list.append(a['element'])
    node_list.append(a['id'])
    weight_list.append(a['weight'])
    atom_names.append(a['atom_name'])

res_dict = {
    'res_name': res_name,
    'nodes': node_list,
    'names': element_list,
    'type':atom_names,
    'weights': weight_list,
    'edges': bond_list}
dict_list.append(res_dict)
    
import pprint
with open('res_dict.txt', 'w', encoding='utf-8') as f:
    # sort_dicts=False keeps keys in the same order as written
    # width=200 ensures lines have room to stretch out before wrapping
    pprint.pprint(dict_list, f, sort_dicts=False, width=300)
