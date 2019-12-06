#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io
import pickle
from functions import identify_N
import os

def qm_region(data,traj, struc_dir,  N_list, total_original_atoms):
	'''
	identifies elements in qm region of a zeolite structure based on Al atoms
	inputs : 
		 data  : data dictionary
		 traj  : name of the traj of the structure
		 N_list: neighbor list of all atoms
		total_original_atoms: number of atoms in original zeolite
	outputs: data dictionary with qm region as part of the data[traj_name]
	'''

	data[traj]['qm_region'] =  []			#add entry for qm_region
	atoms     = io.read(struc_dir+'/'+traj)
	neighbors, neighbors['Si'], neighbors['O'] = {},{'N':[]},{'N':[]}

	'identify Al atoms'
	Al_atoms = [] 	
	for index, atom in enumerate(atoms):
		if atom.symbol == 'Al':
			Al_atoms.append(index)

	'identify neighboring Si/O'
	for Al in Al_atoms:
		data[traj]['qm_region'].append(Al)
		neighbors = identify_N(N_list[Al], N_list, neighbors, Al)
	
	'add Si/O elements to qm region'
	for item in neighbors:
		for i in neighbors[item]['N']:
			if i in data[traj]['qm_region']:
				continue
			else:
				data[traj]['qm_region'].append(i)

	'identify O atoms connected to two Si in the qm region, not accounted for yet'

	for item in N_list:
		if set(N_list[item]) <= set(data[traj]['qm_region']):
			if item in data[traj]['qm_region']:
				continue
			else:
				data[traj]['qm_region'].append(item)

	if len(atoms) > total_original_atoms:
		'adding H/metal atoms in QM region'
		for index in range(total_original_atoms, len(atoms)):
			data[traj]['qm_region'].append(index)
				
	return data
