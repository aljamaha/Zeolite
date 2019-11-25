#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io, Atom
import os
from copy import deepcopy
from molmod import *

'''
Generates unique zeolite structure with 1 or 2 Al substituting Si and enumerate adsorption sites
'''

'Inputs'
zeolite = io.read('CHA.cif')	#Zeolite structure
Al  	= 101			#Si to be replaced by Al

'Inputs (dont change)'
cwd  	= os.getcwd()
struc_dir = cwd+'/structures'	#dir to store structures
index 	= 0			#index of the structure
data 	= {}			#store details of each structure
neighbors = {}			#storing neighbors for Si and O
neighbors['O']  = {'N':[],'NN':[],'NNN':[]}
neighbors['Si'] = {'N':[],'NN':[],'NNN':[]}

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

def add_Pd(zeolite, ads, H_position):
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
	zeolite_copy.append(Atom(ads,(H_position[0], H_position[1], H_position[2]+1)))
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

'substitute Si with Al'
zeolite[Al].symbol = 'Al'

'neighbour list'
zeolite.write('tmp.xyz')
N_list = neighbor_list('tmp.xyz')	#dict of neighbors list
os.system('rm tmp.xyz')

'''Building Si and O N, NN, and NNN'''

'Identify N (Si and O)'
for item in N_list[Al]:
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

'Identify NN O'
for item in neighbors['Si']['N']:
	for neighbor in N_list[item]:
		if neighbor in neighbors['O']['N']:
			continue
		elif neighbor in neighbors['O']['NN']:
			continue
		else:
			neighbors['O']['NN'].append(neighbor)

'Identify NN Si'
for item in neighbors['O']['NN']:
	for neighbor in N_list[item]:
		if neighbor in neighbors['Si']['N']:
			continue
		elif neighbor in neighbors['Si']['NN']:
			continue
		else:
			neighbors['Si']['NN'].append(neighbor)

'Identify NNN O'
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

'Identify NNN Si'
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


#for item in data:
#	print(item, data[item])

'''
To do ...
* print further details beside the name in the .traj
* optimal way of adding H/Pd as an adsorbate
* How would things change if Pd had a +1 oxidation state?
* change adsorbate name from Pd to a variable (to accomodate different oxidation states)
* how many total calculations are needed?
* verify neighbouring list always has atoms with 4 or 2 neighbors
   [check every tested Si has 4 neighbors and every tested O has 2 neighbors]
   [when I add a periodic image, am I repeated myself or identifying new combinations?]
   [or check when I use CHA structure by Jeroen]
* when we have 2 Al, there should be structures with 1 Al, so a M+1 could compensate
* change inputs dictionary
* add_Pd must be updated to reflect adsorption not on H
* adding metal to oxidaiton +2 is not the best way here


Questions:
* Do metals repalce Al?
* How can I verify two structures are not symmetric?
* How to optimize placment of H/Pd atoms?
* Should I add a criteria for distance?
* Andrew Gettson [in Ford] did work on Al distribution on chabasize
* am I missing other Pd oxidation states?

identify_repeat_structures(index):
'''
