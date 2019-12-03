#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io, Atom
import os
from copy import deepcopy
from molmod import *
import pickle
from functions import *
from qm_region import qm_region

'''
Generates unique zeolite structure with 1 or 2 Al substituting Si and enumerate adsorption sites
'''

'Inputs'
zeolite = io.read('CHA-T696.xyz')	#Zeolite structure
Al	= 0

'Inputs (dont change)'
cwd  	= os.getcwd()
struc_dir = cwd+'/structures'	#dir to store structures
index 	= 0			#index of the structure
data 	= {}			#store details of each structure
neighbors = {}			#storing neighbors for Si and O
neighbors['O']  = {'N':[],'NN':[],'NNN':[]}
neighbors['Si'] = {'N':[],'NN':[],'NNN':[]}

'substitute Si with Al'
zeolite[Al].symbol = 'Al'

'neighbour list'
zeolite.write('tmp.xyz')
N_list = neighbor_list('tmp.xyz')	#dict of neighbors list
os.system('rm tmp.xyz')

pickle.dump(N_list, open("save.p", "wb"))

'''Building Si and O N, NN, and NNN'''
neighbors = identify_N(N_list[Al],N_list, neighbors, Al)
neighbors = identify_NN_O(neighbors, N_list)
neighbors = identify_NN_Si(neighbors, N_list)
neighbors = identify_NNN_O(N_list, neighbors)
neighbors = identify_NNN_Si(N_list, neighbors)

'''Writing structures [one Al, two Al [NN and NNN]]'''
'single Al'
index, data = print_structure(zeolite, index, 'N', str(index+1)+'.traj' , struc_dir, data)

'2 Al [NN]'
for item in neighbors['Si']['NN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index, data = print_structure(zeolite_copy, index, 'NN', str(index+1)+'.traj',struc_dir, data)

'2 Al [NNN]'
for item in neighbors['Si']['NNN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index, data = print_structure(zeolite_copy, index, 'NNN', str(index+1)+'.traj',struc_dir, data)

'''Writing structures of H-zeolites'''
zeolite_bare = list(data.keys())	#list of zeolites with Al but no H
index, data  = H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list)

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
							index, data = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'],struc_dir, data)
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
							index, data = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'],struc_dir, data)
							break

'''identify qm region'''
for item in data:
	data = qm_region(data, item, struc_dir, N_list)
	print(data[item])

'''
To do ...
* enough space to accommodate for NO?
* when I add a periodic image, am I repeated myself or identifying new combinations?]
* adding metal to oxidaiton +2 adds only one metal on one of the two Al sites
* am I adding Pd to the optimal site?
* am I missing other Pd oxidation states?

Later ...
* change inputs dictionary
* change adsorbate name from Pd to a variable (to accomodate different oxidation states)

Questions:
* How can I verify two structures are not symmetric?
* Should I add a criteria for distance?
* Andrew Gettson [in Ford] did work on Al distribution on chabasize

*** Other parts of the code: ***
identify_repeat_structures(index)d = []
for i in range(1,index):
	d.append(Al_Al_distance(io.read(struc_dir+'/'+str(i)+'.traj')))

Finds the distance between metal atoms and nearby elements
d = {}
for item in metal_zeolites:
	atoms = io.read(struc_dir+'/'+item)
	atoms.write('tmp.xyz')
	N_list = neighbor_list('tmp.xyz')	#dict of neighbors list
	os.system('rm tmp.xyz')
	d[item] = []
	for N in N_list[metal_zeolites[item]]:
		d[item].append(atoms.get_distance(metal_zeolites[item],N))

Space surrounding metal zeolites

metal_zeolites = {} #store zeolites with metals

for item in data:
	if 'Pd' in data[item].keys():
		atoms = io.read(struc_dir+'/'+item)
		for atom in atoms:
			if atom.symbol == 'Pd':
				metal_zeolites[item] = atom.index
				break
'''
