from ase import io
import pickle, os
from functions import *
from copy import deepcopy

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

def fake_MR(sub_list, N_list):
	'''
	checks that items in the MR do not have more than two neighbors within the MR
	e.g. two 4 MR can make a 6 MR. This should not count as a 6 MR
	Inputs:
		sub_list - list of atoms in the MR
		N_list	 - global neighbor list
	Output:
		True   - not a fake MR
		False  - fake MR that is based on a smaller MR
	'''

	tmp = False
	for k in sub_list:
		n_k = individual_NL(k, N_list)['Si']['N']
		mutual = set(n_k) & set(sub_list)
		if len(mutual) != 2:
			tmp = True

	return tmp

def MR_Al_Al_type(n_MR, Al1, Al2, store, data, traj, N_list, qm_add):
	'''
	Add missing atoms in a n_MR containing two Al atoms
	Inputs:
		n_MR  : type of MR (4/5/6)
		store : list of possible connections [Al1, N, NN, ..] starting from Al1
		Al1   : index of first Al atom
		Al2   : index of second Al atom
		data  : global json database
		traj  : name of traj file
		N_list: global neighbor list
		qm_add: Ture - add to qm, False - don't add it to qm  
	Outputs: data - updated data[traj]['qm_region'] with missing atoms in the n_MR MR
	'''

	MR = ''
	for sub_list in store:
		if MR == '':
			'if it already identified Al-Al pair, no need to go over other MR types'
			sub_list = sub_list[0:n_MR+1] 	#sub_list contains elements from Al1 to n_MR+1
			tmp_Al = False			#check Al1 appear twice in the MR
			Al2_appears = False		#check Al2 appears in the MR
			if sub_list[-1] == Al1:
				'check if last element is returning back to original Al atom, closing the MR'
				for index in sub_list[1:-1]:
					if index == Al1:
						'check Al1 atom doesnt appear twice in the MR'
						tmp_Al = True
					if index == Al2:
						'check Al2 appears in the MR'
						Al2_appears = True

				if tmp_Al == False and Al2_appears == True:
					'Al1 does not appear twice in the 6 MR and Al2 is part of the ring'
					if repeated_MR_item(sub_list) == False:
						'make sure no item in the list is repeated (not a full MR)'
						if fake_MR(sub_list, N_list) == False:
							'check no atom in the MR has more than two neighbors'
							if smaller_mr(sub_list, N_list, data, traj)[0] == False:
								'check MR is not composed of a smaller MR'
								for i in sub_list:
									MR  = True
									if qm_add == True:
										if i not in data[traj]['qm_region']:
											data[traj]['qm_region'].append(i)

	if MR == True:
		data[traj]['Al-Al MR'].append(n_MR)

	return data

def MR_single_Al(n_MR, Al, store, data, traj, N_list, qm_add):
	'''
	Add missing atoms in a n_MR containing Al
	Inputs:
		n_MR  : type of MR (4/5/6)
		store : list of possible connections [Al, N, NN, ..] starting from Al
		Al   : index of first Al atom
		data  : global json database
		traj  : name of traj file
		N_list: global neighbor list
		qm_add: Ture - add to qm, False - don't add it to qm  
	Outputs: data - updated data[traj]['qm_region'] with missing atoms in the n_MR MR
	'''
	sub_sub_list = []
	for sub_list in store:
		sub_list = sub_list[0:n_MR+1] 	#sub_list contains elements from Al to n_MR+1
		tmp_Al = False			#check Al appear twice in the MR
		if sub_list[-1] == Al:
			'check if last element is returning back to original Al atom, closing the MR'
			for index in sub_list[1:-1]:
				if index == Al:
					'check Al atom doesnt appear twice in the MR'
					tmp_Al = True

			if tmp_Al == False:
				'Al does not appear twice in the MR'
				if repeated_MR_item(sub_list) == False:
					'make sure no item in the list is repeated (not a full MR)'
					if fake_MR(sub_list, N_list) == False:
						'check no atom in the MR has more than two neighbors'
						if smaller_mr(sub_list, N_list, data, traj)[0] == False:
							sub_list.sort()
							if sub_list not in sub_sub_list:
								sub_sub_list.append(sub_list)
								data[traj]['Al MR'][n_MR] += 1
							'check MR is not composed of a smaller MR'
							for i in sub_list:
								if qm_add == True:
									if i not in data[traj]['qm_region']:
										data[traj]['qm_region'].append(i)

	return data

def smaller_mr(sub_list, N_list, data, traj):
	'''
	identifies if a mr is composed of a smaller mr by identifying a common neighbor among two atoms that
	is not in the mr
	inputs:
		N_list   - global neighbor list
		data	 - global json data list
		traj	 - name of traj file
		sub_list - list of atoms in the mr
		qn_add   - True - add connecting atom to the qm region
	output:
		True  - a smaller mr constitute this MR
		False - no smaller mr consitute this MR
	'''

	atom_neighbors, skip = {},[]
	tmp = False

	for i in sub_list:
		'nieghbor list of each of the sub_list atoms'
		atom_neighbors[i] = individual_NL(i, N_list)['Si']['N']
	
	'find if two atoms in sub_list have a commond neighbor that is not in the sub_list'
	for i in sub_list:
		skip.append(i)	#to avoid double loop
		for j in sub_list:
			if i != j: #avoid comparing it to itself
				if j not in skip: #avoid double loop
					mutual = set(atom_neighbors[i]) & set(atom_neighbors[j]) #identify mutual elements
					for k in mutual: 
						if k not in sub_list:
							'only concerned with atoms not cnosiderd in the MR'
							tmp = True
							return tmp, data
							#if qm_add == True:
							#	if k not in data[traj]['qm_region']:
							#		data[traj]['qm_region'].append(k)
							#		print(traj, 'two 5 MR to make a 6 MR')
							

	return tmp, data


def update_missing_MR(Al1, N_list, data, traj, n_MR_max, Al_atoms, Al2=''):
	'''
	Identifies missing atoms in 6/8 MR not included in qm region (not NN to Al)
	This is done by building a list [Al_index, N, NN, NNN, ..] and identifying if it completes 6/8 MR
	Inputs:
		Al1 - index of first Al atom
		N_list   - global neighbor list
		data	 - global json data list
		traj	 - name of traj file
		Al2	 - index of second Al atom (if only one Al, then Al2='')
		n_MR_max - max number of elements in the MR (usually stops at 8)
	Outputs:
		updated data[traj][qm_atom] with missing atoms in 6/8 MR now included
	'''

	store = {}
	store[Al1],store[Al2] = [],[]	#list of connections from Al_index until x NN [Al, N, NN, NNN, ..]

	'Al neighbor list'
	for al in Al_atoms:
		Al_NL = individual_NL(al, N_list)['Si']['N']
		for index in Al_NL:
			store[al].append([al, index])

		'Develop neighbor lists based on the last item of each list'
		for n in range(0, n_MR_max-1): 
			copy_store, store[al] = deepcopy(store[al]),[]
			for sub_list in copy_store:
				for index in individual_NL(sub_list[-1], N_list)['Si']['N']:
					'add to the list based on the neighbor last item in the sub_list'
					tmp = deepcopy(sub_list)
					tmp.append(index)
					store[al].append(tmp)
	
	'To be adjusted later'
	MR_types = {}
	MR_types['1Al in 1Al']  = [4,5]
	MR_types['1Al in 2Al']  = [4]
	MR_types['Al-Al'] 	= [4]
	
	'single Al MR'
	for al in Al_atoms:
		if len(Al_atoms) == 1:
			mr_single = [4,5]
		else:
			mr_single = [4]
		
		for mr in mr_single:
			data = MR_single_Al(mr, al ,store[al] , data, traj, N_list, True)
		for mr in range(mr_single[-1],7):
			data = MR_single_Al(mr, al, store[al] , data, traj, N_list, False)

		'Al-Al pair MR'
		if len(Al_atoms) > 1:
			for second_al in Al_atoms:
				if al != second_al:
					for mr in range(4,7):
						'adding missing MR in the Al-Al pair'
						data = MR_Al_Al_type(mr, al, second_al, store[al], data, traj, N_list, True)
					#for mr in range(7,13):
					#	'Identifying MR in Al-Al pairs that are > n_MR_max without adding them to qm region'
					#	data = MR_Al_Al_type(mr, al, second_al, store[al], data, traj, N_list, False)
	
	return data

def qm_region(data, traj, struc_dir,  N_list, total_original_atoms, cutoff, n_MR_max, turn_cutoff='off'):
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
		 turn_cutoff : if on, it adds O atoms <cutoff
		n_MR_max - max number of MR to be tested
	Outputs: 
		data   : data dictionary updated with qm region as part of the data[traj_name]
	'''

	data[traj]['qm_region'],data[traj]['Al-Al MR'],data[traj]['Al MR'] = [],[],{} #initiate entries	
	atoms     = io.read(struc_dir+'/'+traj)		#reads the zeolite structure
	neighbors, neighbors['Si'], neighbors['O'] = {},{'N':[]},{'N':[]}

	'preparing dictionary for single Al MR types'
	for j in range(4,n_MR_max+1):
		data[traj]['Al MR'][j] = 0

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

	'identify remaining atoms in the MR'
	if len(Al_atoms) == 2:
		update_missing_MR(Al_atoms[0], N_list, data, traj, n_MR_max, Al_atoms, Al_atoms[1])
	elif len(Al_atoms) == 1:
		update_missing_MR(Al_atoms[0], N_list, data, traj, n_MR_max, Al_atoms, Al2='')

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
