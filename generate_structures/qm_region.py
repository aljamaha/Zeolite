from ase import io
import pickle, os
from functions import *
from copy import deepcopy

def individual_NL(index, N_list):	
	'creates a neighboring list for atom of interest (needs only atom index and complete NL'

	n_list = {} #local neighbor list specific to Si atoms next to O next to terminal Si
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

def MR_4_AlAl(Al1, Al2, N_list, data, traj):
	'''
	Finds if the Al-Al pair is in a 4 MR or not (assumes Al-Al pair in NN position)
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
	outcome = '' 				#assume at the beginning not in a 4 MR
	Al1_n = individual_NL(Al1, N_list)	#neighbor list of Al1
	Al2_n = individual_NL(Al2, N_list)	#neighbor list of Al2

	mutual_Si = list( set(Al1_n['Si']['N']) & set(Al2_n['Si']['N']) ) #find mutual elemments in neighbor

	if len(mutual_Si) == 2:
		'if two mutual atoms in the two Al neighbor list, it means it is in a 4 MR'
		outcome = '4MR'
		print('4 MR')

	return outcome

def MR_4_Si_Si(Si1, Si2, N_list, qm_region):
	'''
	Finds if the Si-Si (next to each other) pair is in a 4 MR and the remaining two items in 4MR are in qm region
	Inputs:
		Si1    - index of Al atom
		Si2    - index of second Al atom
		N_list - dict of complete neighbor list
		qm_region - indicies of atoms in qm region
	Outputs:
		True  - Si-Si pair are in a 4 MR
		False - Si-Si pair not in a 4 MR
	'''

	store = [] 		#list of connections from Al_index until x NN [Al, N, NN, NNN, ..]
	outcome = False 	#assume they do not make a 4 MR

	'building Si1 neighbor list list [si, index]'
	Si1_NL = individual_NL(Si1, N_list)['Si']['N']
	for index in Si1_NL:
		store.append([Si1, index])
	
	'building Si1 neighbor list list [si, index, N, NN, NNN]'
	for n in range(0, 3): 
		copy_store, store = deepcopy(store),[]
		for sub_list in copy_store:
			for index in individual_NL(sub_list[-1], N_list)['Si']['N']:
				'add to the list based on the neighbor last item in the sub_list'
				tmp = deepcopy(sub_list)
				tmp.append(index)
				store.append(tmp)

	for sub_list in store:
		if sub_list[-1] == Si1 and sub_list[1] == Si2:
			'check the first and last element are the same, and Si2 is next to Si1'
			if sub_list[2] in qm_region and sub_list[3] in qm_region:
				'check the two remaining atoms are in qm region'
				if repeated_MR_item(sub_list) == False:
					outcome = True

	return outcome

def MR_4_Si(Si, N_list, data, traj):
	'''
	Finds if Si is in a 4 MR and the is connected to two atoms in qm region
	Inputs:
		Si     - index of Si atom of interest
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
			'check that both are in qm region'
			if all(x in data[traj]['qm_region'] for x in mutual) == True:
				outcome = True

	return outcome

def MR_6(store, Al_index,data, traj, N_list, Al2):
	'''
	Add missing atoms in 6 MR containing two Al atoms
	Inputs:
		store: list of possible connections [Al, N, NN, ..] that would for a 6 MR
		Al_index: index of first Al atom
		data: global json database
		traj: name of traj file
		N_list: global neighbor list
		Al2: index of second Al atom
	Outputs: data - updated data[traj]['qm_region'] with missing atoms in 6 MR
	'''

	MR = ''
	for sub_list in store:
		sub_list = sub_list[0:7] 	#since we are checking for a 6 MR
		tmp_Al = False			#check Al1 appear twice in the 6 MR	
		Al2_appears = False		#check Al2 appear in the 6 MR
		if sub_list[-1] == Al_index:
			'check if last element is returning back to original Al atom, closing a 6 MR'
			for index in sub_list[1:-1]:
				if index == Al_index:
					'check Al atom doesnt appear twice in a 6 MR (mistaken 6 MR for an extended 4 MR)'
					tmp_Al = True
				if index == Al2:
					'check Al2 appears in the 6 MR'
					Al2_appears = True

			if tmp_Al == False and Al2_appears == True:
				'Al1 does not appear twice in the 6 MR and Al2 is part of the ring'
				if repeated_MR_item(sub_list) == False:
					'make sure no item in the list is repeated (not a full MR)'
					for i in sub_list:
						if i not in data[traj]['qm_region']:
							'add it only if not in the qm region'
							if MR_4_Si(i, N_list, data, traj) == False:
								'check not a member of a 4MR'
								data[traj]['qm_region'].append(i)
								print('6 MR!', i, sub_list)
								MR  = '6MR'

	return data, MR

def MR_8(store, Al_index, Al2, data, traj, N_list):
	'''
	Identifies and adds missing atoms in 8 MR containing two Al atoms	
	Inputs:
		store: list of possible connections [Al, N, NN, ..] that would for a 6 MR
		Al_index: index of first Al atom
		data: global json database
		traj: name of traj file
		N_list: global neighbor list
		Al2: index of second Al atom
	Outputs: data - updated data[traj]['qm_region'] with missing atoms in 8 MR
	'''

	candidates = {}
	for sub_list in store:
		tmp_Al = False		#check Al1 appear twice in the 8 MR	
		Al2_appears = False	#check Al2 appear in the 8 MR
		if sub_list[-1] == Al_index:
			'check if last element is returning back to original Al atom, closing a 8 MR'
			for index in sub_list[1:-1]:
				if index == Al_index:
					'check Al atom doesnt appear twice in a 8 MR (mistaken 8 MR for an extended 4 MR)'
					tmp_Al = True
				if index == Al2:
					'check Al2 appears in the 8 MR'
					Al2_appears = True
			if tmp_Al == False and Al2_appears == True:
				'Al1 does not appear twice in the 8 MR and Al2 is part of the ring'
				if repeated_MR_item(sub_list) == False:
					'make sure no item in the list is repeated (not a full MR)'
					for i in sub_list:
						if i not in data[traj]['qm_region']:
							'add it only if not in the qm region'
							if MR_4_Si(i, N_list, data, traj) == False:
								'check not a member of a 4MR'
								candidates[i] = sub_list

	'check for candidates to be added, they have two connections to the qm region'
	'eliminates dangling Si'
	'first test ...'
	qm_candidates = list(candidates.keys()) + data[traj]['qm_region']
	candidates2 = []
	for Si in list(candidates.keys()):
		Si_n = individual_NL(Si, N_list)['Si']['N']
		n_in_qm = 0
		for neighbor in Si_n:
			if neighbor in qm_candidates:
				n_in_qm += 1
		if n_in_qm > 1:
			if Si not in data[traj]['qm_region']:
				candidates2.append(Si)

	'second test ..'	
	qm_candidates = candidates2 + data[traj]['qm_region']
	candidates3 = []
	for Si in candidates2:
		Si_n = individual_NL(Si, N_list)['Si']['N']
		n_in_qm = 0
		for neighbor in Si_n:
			if neighbor in qm_candidates:
				n_in_qm += 1
		if n_in_qm > 1:
			if Si not in data[traj]['qm_region']:
				candidates3.append(Si)

	'third test ..'
	candidates4 = []
	for Si in candidates3:
		if all(x in data[traj]['qm_region'] + candidates3 for x in candidates[Si]) == True:
			candidates4.append(Si)

	'fourth test ..'
	no_add = []
	for item1 in candidates4: 
		for item2 in candidates4:
			if item1 != item2:
				results = MR_4_Si_Si(item1, item2, N_list, candidates4+data[traj]['qm_region'])
				if results == True:
					if item1 not in no_add:
						no_add.append(item1)
					if item2 not in no_add:
						no_add.append(item2)

	'fifth test ...'
	for k in candidates4:
		if k not in no_add:
			data[traj]['qm_region'].append(k)
			print('8 MR!', k)
			
	return data

def repeated_MR_item(MR_list):
	'''
	Checks if an item in the MR is repeated
	Inputs: 
		MR_list - list of the item in the MR
		
	Output:
		True  - there is a repeated item in the list
		False - there is no repeated item in the list
	'''

	outcome = False
	MR_list = MR_list[1:] #since first and last element should be the same

	a = set([x for x in MR_list if MR_list.count(x) > 1])
	if len(a) > 0:
		outcome = True

	return outcome

def update_missing_MR(Al_index, N_list, data, traj, Al2, n_MR_max):
	'''
	Identifies missing atoms in 6/8 MR not included in qm region (not NN to Al)
	This is done by building a list [Al_index, N, NN, NNN, ..] and identifying if it completes 6/8 MR
	Inputs:
		Al_index - index of first Al atom
		N_list   - global neighbor list
		data	 - global json data list
		traj	 - name of traj file
		Al2	 - index of second Al atom
		n_MR_max - max number of elements in the MR (usually stops at 8)
	Outputs:
		updated data[traj][qm_atom] with missing atoms in 6/8 MR now included
	'''

	store = [] 	#list of connections from Al_index until x NN [Al, N, NN, NNN, ..]
	original_qm = deepcopy(data[traj]['qm_region'])	#qm region before any atom is added

	'Al neighbor list list'
	Al_NL = individual_NL(Al_index, N_list)['Si']['N']
	for index in Al_NL:
		store.append([Al_index, index])
	
	'Develop neighbor lists based on the last item of each list'
	for n in range(0, n_MR_max-1): 
		copy_store, store = deepcopy(store),[]
		for sub_list in copy_store:
			for index in individual_NL(sub_list[-1], N_list)['Si']['N']:
				'add to the list based on the neighbor last item in the sub_list'
				tmp = deepcopy(sub_list)
				tmp.append(index)
				store.append(tmp)

	'check Al-Al not in a 4MR'	
	MR = MR_4_AlAl(Al_index, Al2, N_list, data, traj)

	if MR != '4MR':
		'make sure Al-Al pair not in a 4 MR, then check if it is in a 6 MR'
		data, MR = MR_6(store, Al_index,data, traj, N_list, Al2)

		if MR != '6MR': 
			data = MR_8(store,  Al_index, Al2, data, traj, N_list)

	return data
				
def qm_region(data, traj, struc_dir,  N_list, total_original_atoms):
	'''
	Objective:
		- identifies elements in qm region of a zeolite structure based on Al atoms
		- includes Al atoms, neighboring Si to the Al atoms, and oxygens connecting Si to Al
		- it also includes any adsorbates or protons in the structure
		- includes remianing Si in the 6/8 MR [depending on the case, additional atoms might be added]
	Inputs: 
		 data  : input data dictionary [containing information on each traj file calculation]
		 traj  : name of the traj of the structure
		 N_list: neighbor list of all atoms
		 total_original_atoms: number of atoms in original zeolite (without any protons or adsorbates)
	Outputs: 
		data   : data dictionary updated with qm region as part of the data[traj_name]
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
		update_missing_MR(Al_atoms[0], N_list, data, traj, Al_atoms[1], n_MR_max = 8)

	'identify O atoms connected to Si in the qm region (only ones not accounted for yet)'
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
