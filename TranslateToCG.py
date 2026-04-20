#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 03:32:49 2026

@author: sebas_fu
"""

# code to translate atoms ids from an .ipt or .top file,
# using graph theory, to a CG beads dictionary

# note: atom names can be removed form dictionaries (they are not used),
# however, they help keeping track on the atoms and their weight

import networkx as nx
import pprint
import ast # to read .txt as a list of dictionaries 

# variables
input_file='input_file' #( .mol2, .itp or .top)


# ________________________ Definitions ________________________ 

# parsing of input file (either .itp or .top)
def parse_atoms(line):
    if not line.startswith(';'):
        atom=line.split()
        atom_info = {
            'id': int(atom[0]),
            'element': atom[4][0],
            'res_id': int(atom[2]),
            'res_name': atom[3],
            'weight':round(float(atom[7]))}
        return(atom_info)

def bonds_prase(line):
    if not line.startswith(';'):
        bond=line.split()
        bond_info = {
            'a_i': int(bond[0]),
            'a_j': int(bond[1])}
        return(bond_info)

# parsing of .mol2 files
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

# to determine non-overlapping (unique) matches
def find_uniques_matches(matches_list, total_nodes, ignore_total_nodes=True):
    # if matching doesn't contain all atoms of the residue, change ignore_total_nodes=True

    # Create a compatibility graph to determine overlapping
    compatibility_graph = nx.Graph()

    # Add nodes to compatibility graph
    for i in range(len(matches_list)):
        compatibility_graph.add_node(i)

    # Connect two matches if they don't overlap
    for i in range(len(matches_list)):
        nodes_i = set(matches_list[i].keys()) # for each match, a set of its atoms ids

        for j in range(i + 1, len(matches_list)):
            nodes_j = set(matches_list[j].keys()) # for each other match, other set of its atoms ids

            # If they don't share atoms. Connect them.
            # Makes a connection between each match that with different atoms
            if nodes_i.isdisjoint(nodes_j):
                compatibility_graph.add_edge(i, j)

    # Find cliques (overlapping)
    # Clique is group of nodes interconnected to everyone.
    # This means that all matched in cliques have diffetent atoms (thus, no overlapping)
    all_cliques = list(nx.find_cliques(compatibility_graph))

    # Convert indices back to matching data
    solutions = []
    for clique in all_cliques:
        solution_set = [matches_list[i] for i in clique]
        solutions.append(solution_set)

    # filter solutions so n_nodes = n_total_nodes
    if ignore_total_nodes==False:
        final_solutions=[]
        for sol in solutions:
            n_nodes = sum(len(m) for m in sol) # count nodes in solutions
            if n_nodes==total_nodes:
                final_solutions.append(sol)
        solutions=final_solutions
    return solutions

def create_graph(node_list, name_list, edge_list, weight_list):
    graph_ref = nx.Graph() # backbone
    dict_ref=[]
    graph_ref.add_edges_from(edge_list)
    for n, nm, w in zip(node_list, name_list,weight_list):
        dict_ref.append({'id': n, 'name':nm})
        graph_ref.add_edge(n,-w) # adds atom mass as an edge between atom and its mass (we use neg value to avoid atom id matching with weight node)
    return(graph_ref, dict_ref)

# ________________________ open files ________________________ 

# input file
with open(f"{input_file}", 'r', encoding="utf-8") as f:
    lines = f.readlines()

# dictionaries
with open('dictionaries/res_dictionary.txt', 'r', encoding='utf-8') as f:
    ref_dict = ast.literal_eval(f.read())

with open('dictionaries/CG_dictionary.txt', 'r', encoding='utf-8') as f:
    CG_dict = ast.literal_eval(f.read())


# warning list
warnings=[]

# parsing
atom_data=[]
bond_data=[]
section=None

if input_file.endswith('.top') or input_file.endswith('.itp'):
    for line in lines:
        line = line.strip() # eliminate space between lines
        if line: # ignores empty lines
            if line.startswith("[ "):
                section = line
                continue # this ignores section line
            
            # Atoms data prasing
            if section == "[ atoms ]":
                a_info=parse_atoms(line)
                if a_info:
                    atom_data.append(a_info)
            
            if section == "[ bonds ]":
                b_info=bonds_prase(line)
                if b_info:
                    bond_data.append(b_info)

if input_file.endswith('.mol2'): # parsing if mol2
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
else:
    warn=f'File format for {input_file} not recognized'
    print(warn)
    warnings.append(warn)

# ________________________ graph mapping ________________________ 
matches_in_total=[]
for n in range(1,atom_data[-1]['res_id']+1): # for each residue
    
    graph_ref = nx.Graph() # to avoid error if first recidue is not in reference
    graph_res = nx.Graph()

    res_atom_ids=[] # nodes for each residue
    res_bonds=[] # list of bonds inside each residue
    res_dict=[] # dictionaries of atoms inside each residue
    res_weights=[]
    res_names=[]

    # filter atoms for each residue
    for a in atom_data:
        if (a['res_id'] == n and a['element'] != 'H'): # Ignore Hydrogens, so graph comparison is less expensive
            res_atom_ids.append(a['id'])
            res_weights.append(a['weight'])
            res_names.append(a['element'])
            res_dict.append(a)
            
        
        # this is done when residue need hydrogens to differenciate atoms (e.g. Serine):
        if (a['res_id'] == n and a['res_name'].startswith('SER') and a['element'] == 'H'):
            res_atom_ids.append(a['id'])
            res_weights.append(a['weight'])
            res_names.append(a['element'])
            res_dict.append(a)
        
    n_atoms=len(res_atom_ids) # numbr of atoms in residue

    # finds bonds involved with atoms inside each residue
    for b in bond_data:
        if b['a_i'] in res_atom_ids and b['a_j'] in res_atom_ids:
            res_bonds.append((b['a_j'], b['a_i']))
    
    # if there are no atoms found in a residue
    if n_atoms==0:
        warn=f'res id {n} has no atoms'
        print(warn)
        warnings.append(warn)

    # create graph of residue    
    if n_atoms != 0:
        matching=False
        res_name=res_dict[0]['res_name']
        for res in ref_dict:
            ref_name=res['res_name']
            if res_name.startswith(f'{ref_name}'):
                matched_res=ref_name
                graph_res, dict_res = create_graph(res_atom_ids, res_names, res_bonds, res_weights)
                graph_ref, dict_ref=create_graph(res['nodes'], res['names'], res['edges'], res['weights'])
                matching=True
                break  # Stop searching
        if not matching:
            warn=f'Warning: residue name {res_name} for res_id {n} was not found in reference dictionary'
            print(warn)
            warnings.append(warn)

        # SubGraph
        matches_raw=[]
        ISMAGS = nx.isomorphism.ISMAGS(graph_res, graph_ref) # importa qué va primero
        LCIS = list(ISMAGS.largest_common_subgraph(symmetry=True)) # Por alguna razón no funcionada con symmetry = True
        
        # Filter to also match chemical elements
        for d in LCIS:
            ele_ref=[]
            ele_res=[]
            for key, value in d.items():
                ele_ref.append(next((item_ref['name'] for item_ref in dict_ref if item_ref['id'] == value), None))
                ele_res.append(next((item['element'] for item in res_dict if item['id'] == key), None))
            if ele_ref == ele_res:
                matches_raw.append(d)
        
        # searching non-overlapping solutions
        u_matches=find_uniques_matches(matches_raw, len(graph_res.nodes()))

        if len(u_matches) > 1:
            warn=f'Warning: {len(u_matches)} matches found for {res_name} {n}'
            print(warn)
            warnings.append(warn)
        if len(u_matches) == 0:
            warn=f'Warning: no matches found for {res_name} {n}'
            print(warn)
            warnings.append(warn)
        
        # For the terminal residue:
        if (n == atom_data[-1]['res_id']):
            if len(u_matches) > 0:
                u_matches=[u_matches[0]] # selects only one match
                warn=f'Residue {res_name} {n} is terminal res, only first matching is used'
                print(warn)
                warnings.append(warn)
        
        # appends unique solutions found
        for u in u_matches:
            matches_in_total.append({'res_name':matched_res, 'dictionary':u, 'res_id':n})

# _____________________ matching with CG ids _____________________ 
final_dict=[]
id_CG=1 # initial id for the CG particle
for m in matches_in_total:
    for d_CG in CG_dict:
        if m['res_name']==d_CG['res_name']:
            for node,name in zip(d_CG['nodes'],d_CG['names']):
                atom_id_list=[]
                for v, k in zip(m['dictionary'][0].values(),m['dictionary'][0].keys()):
                    for n in node:
                        if n == v:
                            atom_id_list.append(k)
                final_dict.append({'CG_id':id_CG,'atom_ids':atom_id_list,
                                   'CG_type':name,'res_name':d_CG['res_name'], 'res_id':m['res_id']})
                id_CG+=1


# ________________________ Output files ________________________ 

# write warnings file
with open('warnings.txt', 'w') as f:
        for w in warnings:
            f.write(w + '\n')

# write coarse grain list of beads on a .txt file
with open('final_dictionary.txt', 'w', encoding='utf-8') as f:
    # sort_dicts=False keeps keys in the same order as written
    # width=200 ensures one line per bead id
    pprint.pprint(final_dict, f, sort_dicts=False, width=300)

