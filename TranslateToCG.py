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

# variables
input_file='secretasa_aa.mol2'

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

# parsing .mol2 files
weight_dict={'H':1, 'C':12, 'N':14, 'O':16, 'S':32} # needed because mol2 files dont have atoms weight

def parse_atoms_mol2(line):
    atom=line.split()
    atom_info = {
        'id': int(atom[0]),
        'element': atom[5][0],
        'res_id': int(atom[6]),
        'res_name': atom[7] if len(atom) > 6 else 'Unk',
        'weight':weight_dict[atom[5][0]]}
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

# residue dictionary
ref_dict=[
    {'res_name':'ALA','nodes':[1,2,3,4,5],
     'names':['N','C','C','C','O'],
     'weights':[14,12,12,12,16],
     'edges':[(1,2),(2,3),(2,4),(4,5)]},
    {'res_name':'ASP','nodes':[1,2,3,4,5,6,7,8],
     'names':['N','C','C','O','C','C','O','O'],
     'weights':[14,12,12,16,12,12,16,16],
     'edges':[(1,2),(2,3),(2,5),(3,4),(5,6),(6,7),(6,8)]},
    {'res_name':'ASN','nodes':[1,2,3,4,5,6,7,8],
     'names':['N','C','C','C','N','O','C','O'],
     'weights':[14,12,12,12,14,16,12,16],
     'edges':[(1,2),(2,3),(2,7),(3,4),(4,5),(4,6),(7,8)]},
    {'res_name':'ARG','nodes':[1,2,3,4,5,6,7,8,9,10,11],
     'names':['N','C','C','O','C', 'C', 'C','N','C','N','N'],
     'weights':[14,12,12,16,12,12,12,14,12,14,14],
     'edges':[(1,2),(2,3),(2,5),(3,4),(5,6),(6,7),(7,8),(8,9),(9,10),(9,11)]},
    {'res_name':'CYS','nodes':[1,2,3,4,5,6],
     'names':['N','C','C','S','C','O'],
     'weights':[14,12,12,32,12,16],
     'edges':[(1,2),(2,3),(2,5),(3,4),(5,6)]},
    {'res_name':'GLU','nodes':[1,2,3,4,5,6,7,8,9],
     'names':['N','C','C','O','C','C','C','O','O'],
     'weights':[14,12,12,16,12,12,12,16,16],
     'edges':[(1, 2), (2, 3), (2, 5), (3, 4), (5, 6), (6, 7), (7, 8), (7, 9)]},
    {'res_name':'GLN','nodes':[1,2,3,4,5,6,7,8,9],
     'names':['N', 'C', 'C', 'C', 'C', 'N', 'O', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 14 , 16 , 12 , 16],
     'edges':[(1, 2), (2, 3), (2, 8), (3, 4), (4, 5), (5, 6), (5, 7), (8, 9)]},
    {'res_name':'GLY','nodes':[1,2,3,4],
     'names':['N', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 16],
     'edges':[(1, 2), (2, 3), (3, 4)]},
    {'res_name':'HIS','nodes':[1,2,3,4,5,6,7,8,9,10],
     'names':['N', 'C', 'C', 'C', 'C', 'N', 'C', 'N', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 14 , 12 , 14 , 12,  16],
     'edges':[(1, 2), (2, 3), (2, 9), (3, 4), (4, 5), (5, 8), (4, 6), (6, 7), (7, 8), (9, 10)]},
    {'res_name':'LEU','nodes':[1,2,3,4,5,6,7,8],
     'names':['N', 'C', 'C', 'C', 'C', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 12 , 12 , 16],
     'edges':[(1, 2), (2, 3), (2, 7), (3, 4), (4, 5), (4, 6), (7, 8)]},
    {'res_name':'MET','nodes':[1,2,3,4,5,6,7,8],
     'names':['N', 'C', 'C', 'C', 'S', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 32 , 12 , 12 , 16],
     'edges':[(1, 2), (2, 3), (2, 7), (3, 4), (4, 5), (5, 6), (7, 8)]},
    {'res_name':'PHE','nodes':[1,2,3,4,5,6,7,8,9,10,11],
     'names':['N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 12 , 12 , 12, 12 , 12 , 16 ],
     'edges':[(1, 2), (2, 3), (2, 10), (3, 4), (4, 5), (4, 6), (6, 8), (5, 7), (8, 9), (7, 9), (10, 11)]},
    {'res_name':'PRO','nodes':[1,2,3,4,5,6,7],
     'names':['N', 'C', 'C', 'C', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 12 , 16 ],
     'edges':[(1, 2), (1, 5), (2, 3), (2, 6), (3, 4), (4, 5), (6, 7)]},
    {'res_name':'SER','nodes':[1,2,3,4,5,6,7,8,9,10,11], # includes hydrogen
     'names':['N', 'C', 'C', 'O', 'C', 'O', 'H', 'H', 'H', 'H', 'H'],
     'weights':[14, 12 , 12 , 16 , 12 , 16 , 1  , 1  , 1  , 1  , 1 ],
     'edges':[(1, 2), (1, 9), (2, 3), (2, 5), (2, 10), (3, 4), (5, 6), (3, 7), (3, 8), (4, 11)]},
    {'res_name':'THR','nodes':[1,2,3,4,5,6,7],
     'names':['N', 'C', 'C', 'C', 'O', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 16 , 12 , 16],
     'edges':[(1, 2), (2, 3), (2, 6), (3, 4), (3, 5), (6, 7)]},
    {'res_name':'TRP','nodes':[1,2,3,4,5,6,7,8,9,10,11,12,13,14],
     'names':['N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'O'],
     'weights':[14, 12 , 12 , 12 , 12 , 12 , 12, 12 , 14 , 12 , 12 , 12 , 12 , 16 ],
     'edges':[(1, 2), (2, 3), (2, 13), (3, 4), (4, 6), (4, 5), (6, 8), (6, 7), 
              (7, 9), (8, 11), (5, 9), (7, 10), (11, 12), (10, 12), (13, 14)]},
    {'res_name':'ILE','nodes':[1,2,3,4,5,6,7,8],
     'names':['N','C','C','C','C','C','C','O'],
     'weights':[14,12,12,12,12,12,12,16],
     'edges':[(1,2),(2,3),(2,7),(3,4),(3,5),(5,6),(7,8)]},
    {'res_name':'VAL','nodes':[1,2,3,4,5,6,7],
     'names':['N','C','C','O','C','C','C'],
     'weights':[14,12,12,16,12,12,12,12],
     'edges':[(1,2),(2,3),(2,5),(3,4),(5,6),(5,7)]},
    {'res_name':'TYR','nodes':[1,2,3,4,5,6,7,8,9,10,11,12],
     'names':['N','C','C','C','C', 'C', 'C','C','C','O','C','O'],
     'weights':[14,12,12,12,12,12,12,12,12,16,12,16],
     'edges':[(1,2),(2,3),(2,11),(3,4),(4,5),(4,6),(6,8),(5,7),(8,9),(7,9),(9,10),(11,12)]},
    {'res_name':'LYS','nodes':[1,2,3,4,5,6,7,8,9],
     'names':['N','C','C','O','C', 'C', 'C','C','N'],
     'weights':[14,12,12,16,12,12,12,12,14],
     'edges':[(1,2),(2,3),(2,5),(3,4),(5,6),(6,7),(7,8),(8,9)]}
    ]

CG_dict=[
    {'res_name':'ALA',
    'nodes':[[1,2,3,4],[5]],
    'names':['SP2', 'TC3']},
    {'res_name':'ASP',
    'nodes':[[1,2,3,4],[5,6,7,8]],
    'names':['P2', 'TC3']},
    {'res_name':'ASN',
    'nodes':[[1,2,7,8],[3,4,5,6]],
    'names':['P2', 'SP5']},
    {'res_name':'ARG',
    'nodes':[[1,2,3,4],[5,6,7],[8,9,10,11]],
    'names':['P2', 'SC3','SQ3p']},
    {'res_name':'CYS',
    'nodes':[[1,2,5,6],[3,4]],
    'names':['P2', 'TC6']},
    {'res_name':'GLU',
    'nodes':[[1,2,3,4],[5,6,7,8,9]],
    'names':['P2', 'Q5n']},
    {'res_name':'GLN',
    'nodes':[[1,2,8,9],[3,4,5,6,7]],
    'names':['P2', 'P5']},
    {'res_name':'GLY',
    'nodes':[[1,2,3,4]],
    'names':['SP1']},
    {'res_name':'HIS',
    'nodes':[[1,2,9,10],[3,4],[5,8],[6,7]],
    'names':['P2', 'TC4', 'TN6d', 'TN5a']},
    {'res_name':'LEU',
    'nodes':[[1,2,7,8],[3,4,5,6]],
    'names':['P2', 'SC2']},
    {'res_name':'MET',
    'nodes':[[1,2,7,8],[3,4,5,6]],
    'names':['P2', 'C6']},
    {'res_name':'PHE',
    'nodes':[[1,2,10,11],[3,4],[5,7,9],[6,8]],
    'names':['P2', 'SC4', 'TC5', 'TC5']},
    {'res_name':'PRO',
    'nodes':[[1,2,6,7],[3,4,5]],
    'names':['SP2a', 'SC3']},
    {'res_name':'SER',
    'nodes':[[1,2,5,6],[3,4]],
    'names':['P2', 'TP1']},
    {'res_name':'THR',
    'nodes':[[1,2,6,7],[3,4,5]],
    'names':['P2', 'SP1']},
    {'res_name':'TRP',
    'nodes':[[1,2,13,14],[3,4],[5,9],[6,7],[8,11],[10,12]],
    'names':['P2', 'TC4', 'TN6d', 'TC5', 'TC5', 'TC5']},
    {'res_name':'VAL',
    'nodes':[[1,2,3,4],[5,6,7]],
    'names':['SP2', 'SC3']},
    {'res_name':'ILE',
    'nodes':[[1,2,7,8],[3,4,5,6]],
    'names':['SC1', 'P5']},
    {'res_name':'TYR',
    'nodes':[[1,2,11,12],[3,4],[5,7],[6,8],[9,10]],
    'names':['P2', 'TC4', 'TC5', 'TC5', 'TN6']},
    {'res_name':'LYS',
    'nodes':[[1,2,3,4],[5,6,7],[8,9]],
    'names':['P2', 'SC3', 'SQ4p']}
    ]

# reading file
with open(f"{input_file}", 'r', encoding="utf-8") as f:
    lines = f.readlines()

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
    print(f'File format for {input_file} not recognized')

warnings=[]
matches_in_total=[]
for n in range(1,atom_data[-1]['res_id']+1): # for each residue

    graph_res = nx.Graph()
    res_atom_ids=[] # nodes created for each residue
    res_bonds=[] # list of bonds inside each residue
    res_dict=[] # dictionaries of atoms inside each residue
    res_weights=[]
    res_names=[]
    new_hydrogen_id=[]

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
        

    n_atoms=len(res_atom_ids)

    # finds bonds involved with atoms inside each residue
    for b in bond_data:
        if b['a_i'] in res_atom_ids and b['a_j'] in res_atom_ids:
            res_bonds.append((b['a_j'], b['a_i']))
    
    # if there are no atoms found in a residue
    if n_atoms==0:
        print(f'res id {n} has no atoms')

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

    
    
# matching of atom ids with CG ids
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
                                   'CG_name':name,'res_name':d_CG['res_name'], 'res_id':m['res_id']})
                id_CG+=1

# write warnings file
with open('warnings.txt', 'w') as f:
        for w in warnings:
            f.write(w + '\n')

# write coarse grain beads files
with open('CG_dictionary.txt', 'w') as f:
    for d in final_dict:
        f.write(d + '\n')
