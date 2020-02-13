from ase import io
import pickle, os
from functions import *

def individual_NL(index, N_list):	
	'creates a neighboring list for atom of interest (needs only atom index and complete NL'
	n_list = {} #loca neighbor list specific to Si atoms next to O next to terminal Si
	n_list['Si'], n_list['O'] = {},{}
	n_list['Si']['N'], n_list['O']['N'] = [],[]
	n_list['Si']['NN'], n_list['O']['NN'] = [],[]
	n_list['Si']['NNN'], n_list['O']['NNN'] = [],[]
	n_list = identify_N(N_list[index], N_list, n_list, index)
	n_list = identify_NN_O(n_list, N_list )
	n_list = identify_NN_Si(n_list, N_list)
	n_list = identify_NNN_O(N_list, n_list)
	n_list = identify_NNN_Si(N_list, n_list)
	
	return n_list

def MR_4(Al1, Al2, N_list, data, traj):
	'''
	Finds if the Al-Al pair is in a 4 MR or not
	Inputs:
		Al1    - index of Al atom
		Al2    - index of second Al atom
		N_list - dict of complete neighbor list
		data   - json data file for all calculations
		traj   - name of traj file associated with dict
	Outputs:
		4MR - Al-Al pair are in a 4 MR
		''  - Al-Al pair not in a 4 MR
	'''
	outcome = '' 			#assume at the beginning not in a 4 MR
	Al1_n = individual_NL(Al1, N_list)	#neighbor list of Al1
	Al2_n = individual_NL(Al2, N_list)	#neighbor list of Al2

	mutual_Si = list( set(Al1_n['Si']['N']) & set(Al2_n['Si']['N']) ) #find mutual elemments in neighbor

	if len(mutual_Si) == 2:
		'if two mutual atoms in the two Al neighbor list, it means it is in a 4 MR'
		outcome = '4MR'
		print('4 MR')

	return outcome

def MR_4_Si(Si, N_list, data, traj):
	'''finds if Si is in a 4 MR
	Inputs:
		Si    - index of Si atom of interest
		N_list - dict of complete neighbor list
		data   - json data file for all calculations
		traj   - name of traj file associated with dict
	Outputs:
		True: atom is in a 4 MR
		False: atom is not in a 4 MR
	'''

	outcome = False				
	Si_n = individual_NL(Si, N_list)	#neighbor list of Si
	
	for i in Si_n['Si']['NN']:
		'find mutual neighbors between Si and NN'
		test_n = individual_NL(i, N_list)
		mutual = list( set(Si_n['Si']['N']) & set(test_n['Si']['N']) ) #find mutual elemments in neighbor
		if len(mutual) == 2:
			'check that both are in qm region >> to be added'
			if all(x in data[traj]['qm_region'] for x in mutual) == True:
				outcome = True

	return outcome

def middle_connect(Si_number, Al_number, N_list):
	'''
	Finds the middle atom connecting Si_number and Al_number (both are indicies)'
	Input: index of Si atom and Al atom
	Output:
		'' - no single atom between the two atoms
		index of of the middle atom (if exists)
	'''

	output = ''
	Si_list = individual_NL(Si_number, N_list)
	Al_list = individual_NL(Al_number, N_list)
	
	for i in Si_list['Si']['N']:
		if i in Al_list['Si']['N']:
			output = i

	return output

def MR_6_NN(Al1, Al2, N_list, data, traj):	
	'''finds Si atom missing in qm region in 6 MR where Al atoms in NN
	Inputs:
		Al1    - index of Al atom
		Al2    - index of second Al atom
		N_list - dict of complete neighbor list
		data   - json data file for all calculations
		traj   - name of traj file associated with dict
	Outputs:
		adds Si atom (if it is found) to data[traj]['qm_region']
		output - either '6MR' or '' if it is not a 6MR
	'''

	output = ''
	test_list = {}	#dict with entries (atom indexes to be tested) and a list of connecting atoms
	Al1_n = individual_NL(Al1, N_list)	#neighbor list of Al1

	for i in Al1_n['Si']['NN']:
		'start with Al'
		if i not in data[traj]['qm_region']:
			'find NN not in qm region and add it to a list to be checked later'
			test_list[i] = []
			j = middle_connect(Al1, i, N_list)	#Si atom connecting Al and NN Si
			test_list[i] = [Al1, j, i]		#add to the list Al, NN Si and connecting Si

		for element in test_list:
			'check if the Si elment has Al2 as NN. If so, add connecting Si to test list'
			Si = middle_connect(element, Al2, N_list)
			tmp = False
			n   = individual_NL(element, N_list)
			for i in n['Si']['N']:
				if i in data[traj]['qm_region']:
					if Si == i:
						test_list[element].append(i)
						tmp = True
			if tmp == True:
				'if it found NN Al, add Al1 to the test list'
				test_list[element].append(Al2)

		for element in test_list:
			if len(test_list[element]) == 5:
				'it has 5 elements connecting (two Al, 3 Si)'
				'no need to consider Si between two Al because we know it exists (structure is NN)'
				if test_list not in data[traj]['qm_region']: 
					a = MR_4_Si(element, N_list, data, traj)
					if a == False:
						'check it is not in a 4 MR'
						print('6 MR NN', element)
						data[traj]['qm_region'].append(element)
						output = '6MR'

	return data, output

def connecting(data, N_list, traj, atoms, Al1, Al2):
	'''
	Identifies missing atoms connecting 6/8 MR and adds them to QM region
	Inputs:		
		Al1    - index of Al atom
		Al2    - index of second Al atom
		N_list - dict of complete neighbor list
		traj   - name of traj file associated with dict
	Output:
		data - updated with missing qm atoms in 6/8 MR
	'''

	Al_list1 = individual_NL(Al1, N_list)
	Al_list2 = individual_NL(Al2, N_list)

	MR = ''
	MR = MR_4(Al1, Al2, N_list, data, traj)
	
	if MR != '4MR':
		'make sure Al-Al pair not in a 4 MR'
	
		'6 MR'
		data, MR =  MR_6_NN(Al1, Al2, N_list, data, traj)

		'8 MR'
		if MR != '6MR':
			'if it is not a 4/6 MR'
			if data[traj]['N'] == 'NNN':
				'if the NN Si (not in qm region) is shared'
				for a1 in Al_list1['Si']['NN']:
					a1_NL =  individual_NL(a1, N_list)
					for a2 in Al_list2['Si']['NN']:
						if a2 in a1_NL['Si']['N']:
							if a1 not in data[traj]['qm_region'] and a2 not in data[traj]['qm_region']:
								MR = MR_4_Si(a1, N_list, data, traj)
								if MR == False:
									if a1 not in data[traj]['qm_region']:
										data[traj]['qm_region'].append(a1)
									if a2 not in data[traj]['qm_region']:												data[traj]['qm_region'].append(a2)
									print('8 MR NNN')	
								else:
									print('8 is in 4 MR')
		

	return data


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
		
		'6 MR (NN): Do two Al appear twice in NN list? Is it not too close to both Al atoms?'
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
				print('6 MR items added')
				data[traj]['qm_region'].append(item)

		'8 MR (NNN): Do two Al appear in NN and NNN (one in each)'
		if data[traj]['N'] == 'NNN':
			'only use it on NNN'
			Al = {}
			Al['NN'],Al['NNN'] = [],[]
			for n in n_list['Si']['NN']:
				if atoms[n].symbol == 'Al':
					if n not in Al['NN']:
						Al['NN'].append(n)
			for n in n_list['Si']['NNN']:	
				if atoms[n].symbol == 'Al':
					if n not in Al['NNN']:
						Al['NNN'].append(n)
			pass_distance = True
			if len(Al['NN']) == 1 and len(Al['NNN']) == 1:
				for a in Al_neighbor:
					d = atoms.get_distance(a,item)	
					if d < 5.55:
						pass_distance = False
				if pass_distance == True:	
					print('add to 8 MR NNN')
					data[traj]['qm_region'].append(item)
			
		
		'8 MR (NN): Do two Al appear twice in NN list? Is it not too close to both Al atoms?'
		if data[traj]['N'] == 'NN':
			Al_NNN = []
			for n in n_list['Si']['NNN']:
				if atoms[n].symbol == 'Al':
					if n not in Al_NNN:
						Al_NNN.append(n)	
			pass_distance = True
			if len(Al_NNN) == 2:
				for a in Al_neighbor:
					d = atoms.get_distance(a,item)	
					if d < 5.5:
						pass_distance = False
				if pass_distance == True:	
					print('wow!')
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

	'identify remaining atoms in the 6/8 MR'
	if len(Al_atoms) == 2:
		data = connecting(data, N_list, traj, atoms, Al_atoms[0], Al_atoms[1])

	'identify O atoms connected to two Si in the qm region (only ones not accounted for yet)'
	for item in N_list:
		if set(N_list[item]) <= set(data[traj]['qm_region']):
			if item in data[traj]['qm_region']:
				continue
			else:
				data[traj]['qm_region'].append(item)

	'adding H/metal atoms in QM region'
	if len(atoms) > total_original_atoms:
		for index in range(total_original_atoms, len(atoms)):
			data[traj]['qm_region'].append(index)
				
	return data
