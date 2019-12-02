#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io, Atom
import os
from copy import deepcopy
from molmod import *
import pickle

def neighbor_list(xyz_file):
	'''
	Inputs:  xyz coordinates
	Outputs: dictionary of the neighbor list
	'''
	mol = Molecule.from_file(xyz_file)
	mol.set_default_masses()
	mol.set_default_graph()
	return mol.graph.neighbors

def count_elements(atoms):
	'''
	Inputs:  ase atoms object
	Outputs: dictionary of the number of atoms of each element in a structure
	'''
	list_of_symbols, structure_data  = [], {}
	for item in list(atoms.symbols):
		if item not in list_of_symbols:
			list_of_symbols.append(item)
	for item in list_of_symbols:
		structure_data[item] = list(atoms.symbols).count(item)
	return structure_data

def print_structure(atoms, index, N, reference):
	'''
	Inputs:
		atoms: ase atoms object
		index: the index of the previous structure
		N    : N, NN, or NNN
		reference: original structure without modification
	Outputs:
		prints the traj file in struc_dir folder
		returns the index of the traj file in struc_dir
		appends details of the structure to data dictionary
	'''
	index += 1
	os.chdir(struc_dir), atoms.write(str(index)+'.traj')
	data[str(index)+'.traj'] = count_elements(atoms)
	data[str(index)+'.traj']['N'] = N
	data[str(index)+'.traj']['reference'] = reference
	if 'H' not in list(data[str(index)+'.traj']):
		data[str(index)+'.traj']['H'] = 0
	data[str(index)+'.traj']['oxidation']= int(data[str(index)+'.traj']['Al']) - int(data[str(index)+'.traj']['H'])

	return index

def add_H(zeolite, O1, H_num, O2=0):
	'''
	Inputs:
		zeolite: structure of the zeolite
		O      : index of first O atom neighboring Al
		H_num  : number of H atoms in the structure
		O2     : index of second O atom in the second Al
	Outputs:
		zeolite structure with added H atom(s)
	'''
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy.append(Atom('H',(zeolite[O1].position[0], zeolite[O1].position[1], zeolite[O1].position[2]+1)))
	if H_num == 2:
		zeolite_copy.append(Atom('H',(zeolite[O2].position[0],zeolite[O2].position[1],zeolite[O2].position[2]+1)))
	return zeolite_copy

def add_Pd(zeolite, ads, position):
	'''
	Inputs:
		zeolite   : structure of zeolite
		ads       : metal structure to be added
		H_position: xyz coordinates of
	Outputs:
		structure of zeolite with metal adsorbed
	'''
	ads = 'Pd'
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy.append(Atom(ads,(position[0], position[1], position[2]+2)))
	#zeolite_copy.append(Atom(ads,(4.4,5,10)))
	return zeolite_copy

def O_neighbor_indicies(atoms):
	'''
	Inputs:  ase atoms object
	Outputs: dictionary of O atoms neighboring Al[0] and Al[1]
	'''
	Al_index, O_index = [],{}
	for atom in atoms:
		'identify Al atoms'
		if atom.symbol == 'Al':
			Al_index.append(atom.index)
	for i, Al in enumerate(Al_index):
		'identify neighboring O'
		O_index[i] = []
		for O in N_list[Al]:
			O_index[i].append(O)
	return O_index

def identify_repeat_structures(index):
	'''
	Inputs: index
	Output: list of repeated structures
	'''
	output = []
	for i in range(1, index):
		print(i)
		for j in range(1, index):
				#if i != j:
				atoms1 = io.read(str(i)+'.traj')
				atoms2 = io.read(str(j)+'.traj')
				if atoms1 == atoms2:
					print('Same structure! ', str(i)+'\t'+str(j))
					output.append([i,j])
	return output

def Al_Al_distance(atoms):
	'''
	Inputs: atoms object
	Output: Al-Al pair distance
	'''
	n_Al, Al_index = 0,[] #number of Al atoms, index of Al atoms
	for atom in atoms:
		if atom.symbol == 'Al':
			n_Al += 1
			Al_index.append(atom.index)
	if n_Al == 2:
		distance = atoms.get_distance(Al_index[0],Al_index[1])
	else:
		distance = 0

	return distance

def identify_N(N_list_Al, N_list, neighbors, Al):
	'''
	Identify Al neighboring Si and O
	Inputs :
		N_list[Al]: neighbor list of Al
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
		Al 	  : element number of first Al
	Outputs:
		neighbors: Si and O neighboring Al
	'''
	for item in N_list_Al:
		'loop over O attached to Si'
		neighbors['O']['N'].append(item)
		for neighbor in N_list[item]:
			'loop over neighbors of Si attached to neighbor'
			if neighbor == Al:
				continue
			elif neighbor in neighbors['Si']['N']:
				continue
			else:
				neighbors['Si']['N'].append(neighbor)
	return neighbors

def identify_NN_O(neighbors, N_list ):
	'''
	identifies Al next neighbor O
	inputs :
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
	outputs:
		neighbors: O next neighboring Al
	'''
	for item in neighbors['Si']['N']:
		for neighbor in N_list[item]:
			if neighbor in neighbors['O']['N']:
				continue
			elif neighbor in neighbors['O']['NN']:
				continue
			else:
				neighbors['O']['NN'].append(neighbor)
	return neighbors

def identify_NN_Si(neighbors, N_list):
	'''
	identifies Al next neighbor Si
	inputs :
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
	outputs:
		neighbors: Si next neighboring Al
	'''
	for item in neighbors['O']['NN']:
		for neighbor in N_list[item]:
			if neighbor in neighbors['Si']['N']:
				continue
			elif neighbor in neighbors['Si']['NN']:
				continue
			else:
				neighbors['Si']['NN'].append(neighbor)
	return neighbors

def identify_NNN_O(N_list, neighbors):
	'''
	identifies Al next next neighbor O
	inputs :
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
	outputs:
		neighbors: O next next neighboring Al
	'''
	for item in neighbors['Si']['NN']:
		for neighbor in N_list[item]:
			if neighbor in neighbors['O']['N']:
				continue
			elif neighbor in neighbors['O']['NN']:
				continue
			elif neighbor in neighbors['O']['NNN']:
				continue
			else:
				neighbors['O']['NNN'].append(neighbor)
	return neighbors

def identify_NNN_Si(N_list, neighbors):
	'''
	identifies Al next next neighbor O
	inputs :
		N_list    : dictionary of the neighboring list of all atoms
		neighbors : dictionary of neighbors (Si/O)[N,NN,NNN]
	outputs:
		neighbors: Si next next neighboring Al
	'''
	for item in neighbors['O']['NNN']:
		for neighbor in N_list[item]:
			if neighbor in neighbors['Si']['N']:
				continue
			if neighbor in neighbors['Si']['NN']:
				continue
			elif neighbor in neighbors['Si']['NNN']:
				continue
			else:
				neighbors['Si']['NNN'].append(neighbor)
	return neighbors

'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

'''Writing structures [one Al, two Al [NN and NNN]]'''
'single Al'
index = print_structure(zeolite, index, N='N', reference=str(index+1)+'.traj')


'2 Al [NN]'
for item in neighbors['Si']['NN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index = print_structure(zeolite_copy, index, N='NN', reference=str(index+1)+'.traj')

'2 Al [NNN]'
for item in neighbors['Si']['NNN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index = print_structure(zeolite_copy, index, N='NNN', reference=str(index+1)+'.traj')
'''Writing structures of H-zeolites'''
zeolite_bare = list(data.keys())	#list of zeolites with Al but no H

for item in zeolite_bare:
	atoms = io.read(struc_dir+'/'+data[item]['reference'])
	if data[item]['Al'] == 1:
		for O in neighbors['O']['N']:
			zeolite_copy = add_H(atoms, O, 1)
			index = print_structure(zeolite_copy, index, N='N', reference=item)
	else:
		O_index = O_neighbor_indicies(atoms)
		for i in O_index[0]:
			for j in O_index[1]:
				zeolite_copy = add_H(atoms, i, 2, j)
				index = print_structure(zeolite_copy, index, N=data[item]['N'], reference=item)

'''Writing structures of metal modified zeolites'''
no_metal_zeolite = list(data) #List of structures with no introduced metal [includes ones with H]

#### this needs to change####
inputs = {'PdO': [-2, 0, 2], 'Pd2O': [-2, 2, 6], 'PdO2': [-4, -2, 0], 'Pd2O2': [-4, 0, 4], 'Pd': [0, 2, 4]}
ox1, ox2 = 0,0
for structure in no_metal_zeolite:
	if data[structure]['oxidation'] == 0:
		'oxidation state of zero'
		for comp in inputs:
			for ox in inputs[comp]:
				if ox == 0:
					atoms = io.read(struc_dir+'/'+structure)
					for atom in atoms:
						if atom.symbol == 'H':
							zeolite_copy = add_Pd(atoms, comp, atom.position)
							index = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'])
	elif data[structure]['oxidation'] == 1:
		'to be added of +1'
		continue
	elif data[structure]['oxidation'] == 2:
		'oxidation state of +2'
		for comp in inputs:
			for ox in inputs[comp]:
				if ox == 2:
					atoms = io.read(struc_dir+'/'+structure)
					for atom in atoms:
						if atom.symbol == 'Al':
							
							zeolite_copy = add_Pd(atoms, comp, atom.position)
							index = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'])
							break


'''Space surrounding metal zeolites'''

metal_zeolites = {} #store zeolites with metals

for item in data:
	if 'Pd' in data[item].keys():
		atoms = io.read(struc_dir+'/'+item)
		for atom in atoms:
			if atom.symbol == 'Pd':
				metal_zeolites[item] = atom.index
				break

