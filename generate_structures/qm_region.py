#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io
import pickle
from functions import *
import os

def terminal_atoms(data, N_list, traj, atoms):
	'''
	identify remaining atoms in the 6/8 MR
	inputs:
		data
		N_list (neighbor list for all atoms)
		traj  : name of the traj of the structure
		atoms : ase objects atoms
	outputs: 
		updated data with remaining members of 6/8 MR (if available)
	'''
	terminal = {}
	terminal['Si']         = [] #terminal Si in QM region
	terminal['O-Si']       = [] #oxygen connected to terminal Si in QM regiona
	terminal['Si-O-Si(t)'] = [] #Si connected to oxygen connected to terminal Si

	'identify Al atoms'
	Al_atoms = []
	for index, atom in enumerate(atoms):
		if atom.symbol == 'Al':
			Al_atoms.append(index)

	for item in data[traj]['qm_region']:
		'identify terminal Si (connected to only one O in qm region)'
		if atoms[item].symbol == 'Si':
			O_connections = 0
			for n in N_list[item]:
				if n in data[traj]['qm_region']:
					O_connections += 1
			if O_connections == 1:
				terminal['Si'].append(item)

	for item in terminal['Si']:
		'identify O next to terminal Si'
		for n in N_list[item]:
			if n not in terminal['O-Si']:
				terminal['O-Si'].append(n)

	for item in terminal['O-Si']:
		'identify Si connected to oxygen connected to terminal Si'
		for n in N_list[item]:
			if n not in terminal['Si-O-Si(t)']:
				terminal['Si-O-Si(t)'].append(n)

	for item in terminal['Si-O-Si(t)']:
		'Si atom has to be connected to the 6/8 MR to be included in the qm region'	
		n_list = {} #loca neighbor list specific to Si atoms next to O next to terminal Si
		n_list['Si'], n_list['O'] = {},{}
		n_list['Si']['N'], n_list['O']['N'] = [],[]
		n_list['Si']['NN'], n_list['O']['NN'] = [],[]
		n_list['Si']['NNN'], n_list['O']['NNN'] = [],[]
		n_list = identify_N(N_list[item], N_list, n_list, item)
		n_list = identify_NN_O(n_list, N_list )
		n_list = identify_NN_Si(n_list, N_list)
		n_list = identify_NNN_O(N_list, n_list)
		n_list = identify_NNN_Si(N_list, n_list)
		
		'6 MR (NN): Does Al appear twice in NN list? Is it not too close to both Al atoms?'
		Al_neighbor = []
		for n in n_list['Si']['NN']:
			if atoms[n].symbol == 'Al':
				if n not in Al_neighbor:
					Al_neighbor.append(n)
		pass_distance = True
		if len(Al_neighbor) >1:
			'if it passes first test of NN from both Al, check distance from Al atoms is not too short'
			for a in Al_neighbor:
				d = atoms.get_distance(a,item)	
				if d < 5:
					pass_distance = False
			if pass_distance == True:	
				print('YESSS')		
				data[traj]['qm_region'].append(item)
			

	'''	
	for entry in terminal['Si']:

		n_list = {} #neighbor list
		n_list['Si'], n_list['O'] = {},{}
		n_list['Si']['N'], n_list['O']['N'] = [],[]
		n_list['Si']['NN'], n_list['O']['NN'] = [],[]
		n_list['Si']['NNN'], n_list['O']['NNN'] = [],[]
		n_list = identify_N(N_list[entry], N_list, n_list, entry)
		n_list = identify_NN_O(n_list, N_list )
		n_list = identify_NN_Si(n_list, N_list)
		n_list = identify_NNN_O(N_list, n_list)
		n_list = identify_NNN_Si(N_list, n_list)
		
		n_Al = 0
		for j in n_list['Si']['NN']:
			if atoms[j].symbol == 'Al': 
				n_Al += 1
		
		if traj == "3.traj":	
			print('position: ', atoms[entry].position)
			print(entry, n_Al)
			print(n_list)

		if n_Al > 1:
			print('add to qm region')
			data[traj]['qm_region'].append(entry)

	'''
	'''
	for item in terminal['Si']:
		'identify O next to terminal Si'
		for n in N_list[item]:
			if n not in terminal['O-Si']:
				terminal['O-Si'].append(n)
	for item in terminal['O-Si']:
		'identify Si connected to oxygen connected to terminal Si'
		for n in N_list[item]:
			terminal['Si-O-Si(t)'].append(n)

	for item in terminal['Si-O-Si(t)']:
		'if Si occured more than once, it is added '
		n_count = terminal['Si-O-Si(t)'].count(item)
		if n_count > 1:
			if item not in data[traj]['qm_region']:
				data[traj]['qm_region'].append(item)
	'''
	return data

def Al_plane(atoms):
	'''
	Objective: find the axis (plane where the two Al atoms are at
	Inputs: ase aotms object
	plane: plane (xz, xy, yz)
	'''

	'identify Al atoms'
	Al_atoms = []
	for index, atom in enumerate(atoms):
		if atom.symbol == 'Al':
			Al_atoms.append(index)
	
	if len(Al_atoms) == 1:
		plane = ''

	else:
		x = atoms[Al_atoms[0]].position[0] - atoms[Al_atoms[1]].position[0]
		y = atoms[Al_atoms[0]].position[1] - atoms[Al_atoms[1]].position[1]
		z = atoms[Al_atoms[0]].position[2] - atoms[Al_atoms[1]].position[2]

		if x < z and y < z :
			plane =  'xy'
		elif x < y and z < y :
			plane = 'xz'
		elif z < x and y < x:
			plane =  'yz'

	return plane

def qm_region(data, traj, struc_dir,  N_list, total_original_atoms):
	'''
	Objective:
		- identifies elements in qm region of a zeolite structure based on Al atoms
		- includes Al atoms, neighboring Si to the Al atoms, and oxygens connecting Si to Al
		- it also includes any adsorbates or protons in the structure
		- (in progress) includes remianing Si in the 6/8 MR
	Inputs: 
		 data  : input data dictionary [containing information on each traj file calculation]
		 traj  : name of the traj of the structure
		 N_list: neighbor list of all atoms
		 total_original_atoms: number of atoms in original zeolite (without any protons or adsorbates)
	Outputs: data dictionary with qm region as part of the data[traj_name]
	'''

	data[traj]['qm_region'] =  []			#initiate entry for qm_region
	atoms     = io.read(struc_dir+'/'+traj)		#reads the zeolite structure
	neighbors, neighbors['Si'], neighbors['O'] = {},{'N':[]},{'N':[]}

	'identify Al atoms'
	Al_atoms = [] 	
	for index, atom in enumerate(atoms):
		if atom.symbol == 'Al':
			Al_atoms.append(index)

	'identify neighboring Si/O to Al atom(s)'
	for Al in Al_atoms:
		data[traj]['qm_region'].append(Al)
		neighbors = identify_N(N_list[Al], N_list, neighbors, Al)
	
	'add Si/O neighboring Al atom(s) to qm region'
	for item in neighbors:
		for i in neighbors[item]['N']:
			if i in data[traj]['qm_region']:
				continue
			else:
				data[traj]['qm_region'].append(i)

	'identify O atoms connected to two Si in the qm region (only ones not accounted for yet)'
	for item in N_list:
		if set(N_list[item]) <= set(data[traj]['qm_region']):
			if item in data[traj]['qm_region']:
				continue
			else:
				print('not sure here', item)
				data[traj]['qm_region'].append(item)

	'identify remaining atoms in the 6/8 MR'
	data = terminal_atoms(data, N_list, traj, atoms)

	'adding H/metal atoms in QM region'
	if len(atoms) > total_original_atoms:
		for index in range(total_original_atoms, len(atoms)):
			data[traj]['qm_region'].append(index)
				
	return data
