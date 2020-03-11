import os, json
from functions import Al_Al_distance
from ase import io
import numpy as np
from copy import deepcopy

def check(item_data, candidate_data, entry, length=False):
	'''
	checks if item_data and candidate_data have the same values of an entry
	If there is a comparison between the length of two entries, length should be "True"
	returns 'pass' if they match, and 'fail' if they dont match
	'''
	tmp = 'fail'
	if length == True:
		if len(item_data[entry]) == len(candidate_data[entry]):
			tmp = 'pass'
	elif item_data[entry] == candidate_data[entry]:
		tmp = 'pass'

	return tmp

def group_check(duplicates, item, candidate, group):
	tmp = 'fail'
	for key in duplicates:
		if item in duplicates[key] and candidate in duplicates[key]:
			tmp = 'pass'
			break
		elif item in duplicates[key] and candidate not in duplicates[key]:
			duplicates[key].append(candidate)
			tmp = 'pass'	
			break
		elif item not in duplicates[key] and candidate in duplicates[key]:
			duplicates[key].append(candidate)
			tmp = 'pass'	
			break
	
	if tmp == 'fail':
		group +=1
		duplicates['group'+str(group)] = [item, candidate]

	return duplicates, group

def compare_Al_distance(item, candidate, struc_dir, cutoff=0.7):	
	
	tmp = 'fail'
	item_atoms      = io.read(struc_dir+'/'+item)
	candidate_atoms = io.read(struc_dir+'/'+candidate)
	item_Al = Al_Al_distance(item_atoms)
	candidate_Al = Al_Al_distance(candidate_atoms)

	d = np.abs(item_Al - candidate_Al)
	if d<cutoff:
		tmp = 'pass'
	return tmp, d

def remove_duplicates(data, struc_dir):
	'''
	Identify structures with duplicates, and then deleting the duplicates from data dictionay
	Inputs:
		data dictionary containing information about all structures
	Output:
	'''
	duplicates, group = {}, 0
	duplicates['group1'] = {}

	for item in data:
		'aggregate duplicates in the same group under duplicates dictionary'
		for candidate in data:
			if item != candidate:
				if check(data[item], data[candidate], 'Al') == 'pass':
					if check(data[item], data[candidate], 'N') == 'pass':
						if check(data[item], data[candidate], 'qm_region', length=True) == 'pass':
							if compare_Al_distance(item, candidate, struc_dir)[0] == 'pass':
								duplicates, group =  group_check(duplicates, item, candidate, group)

	for j in data:
		'for those who do not have duplicates, make them in the same group'
		tmp = 'fail'
		for item in list(duplicates.values()):
			for i in item:
				if i == j:
					tmp = 'pass'

		if tmp == 'fail':
			group += 1
			duplicates['group'+str(group)] = [j]

	'add unique structures to the non_duplicates list'
	non_duplicates = []
	for item in duplicates:
		non_duplicates.append(duplicates[item][0])

	
	'deleting detected duplicate'
	data_copy = deepcopy(data)
	for item in data_copy:
		if item not in non_duplicates:
			del data[item]

	print('Detecting duplicates ...')
	for item in duplicates:
		print(item, duplicates[item])

	return data


