from ase import io
import os, json
from copy import deepcopy
from functions import *

'Inputs'
cutoff = 7

'Update full atoms based on latest run'
os.system('cp ~/qchem_scripts/qchem-to-ase-all-atoms.py .; python qchem-to-ase-all-atoms.py')
atoms = io.read('full-atoms.xyz')

'Al neighbors'
Al_N = [] #Si neighbors of Al atoms
NL   = neighbor_list('full-atoms.xyz')
for atom in atoms:
	if atom.symbol == 'Al':
		new_neighbors = individual_NL(atom.index, NL)['Si']['N']
		for item in new_neighbors:
			Al_N.append(item)

'load data'
with open('dir_data.json','r') as read:
	data = json.load(read)

'Identify N atom'
for atom in atoms:
	if atom.symbol == 'N':
		N_index = atom.index
		break

'delete atoms outside of the cutoff from qm region'
qm_region = deepcopy(data['qm_region'])
for atom in qm_region:
	if atoms.get_distance(atoms[atom].index, N_index) > cutoff and atoms[atom].symbol == 'Si':
		if atom not in Al_N:
			data['qm_region'].remove(atom)

'Add Si atoms within cutoff to qm region'
for atom in atoms:
	if atoms.get_distance(atom.index, N_index) < cutoff:
		if atom.index not in data['qm_region'] and atom.symbol == 'Si':
			data['qm_region'].append(atom.index)

'Add missing oxygens connecting two Si atoms in Qm region'
for atom in atoms:
	if atom.symbol == 'O':
		if atom.index not in data['qm_region']:
			tmp = True
			for n in NL[atom.index]:
				if n not in data['qm_region']:
					tmp = False
			if tmp == True:
				data['qm_region'].append(atom.index)

'Remove dangling oxygen atoms'
qm_region = deepcopy(data['qm_region'])
for item in qm_region:
	if atoms[item].symbol == 'O':
		for n in NL[item]:
			if n not in data['qm_region']:
				data['qm_region'].remove(item)
				break

with open('dir_data.json', 'w') as write:
	json.dump(data, write, indent=4)
