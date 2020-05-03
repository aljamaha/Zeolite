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
			#print('8 MR!', k)

	return data
