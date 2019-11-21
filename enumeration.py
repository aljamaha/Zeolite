from ase import io,Atoms,Atom
import os
from copy import deepcopy
from molmod import *
from ase.build import add_adsorbate

'Inputs'
zeolite = io.read('CHA.cif')	#Zeolite structure
Al  	= 101			#Si to be replaced by Al

'Inputs dont change'
cwd  	= os.getcwd()
index 	= 0	#index of the structure
neighbors = {}	#storing neighbors for Si and O
neighbors['O']  = {'N':[],'NN':[],'NNN':[]}
neighbors['Si'] = {'N':[],'NN':[],'NNN':[]}
struc_dir = cwd+'/structures'	#dir to store structures

def neighbor_list(xyz_file):
	'returns dict of the neighbor list of all atoms in an xyz file'
	mol = Molecule.from_file(xyz_file)
	mol.set_default_masses()
	mol.set_default_graph()	
	os.system('rm '+xyz_file)
	return mol.graph.neighbors

def print_structure(atoms, index):
	'add traj file in structures directory'		
	index += 1
	os.chdir(struc_dir)
	atoms.write(str(index)+'.traj')
	return index

'substitute Si with Al'
zeolite[Al].symbol = 'Al'

'neighbour list'
zeolite.write('a.xyz')
N_list = neighbor_list('a.xyz')	#dict of neighbors list

'''
'Identify and write H atoms next to Al'
for item in N_list[Al]:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy.append(Atom('H',(zeolite[item].position[0], zeolite[item].position[1], zeolite[item].position[2]+1)))
	index = print_structure(zeolite_copy, index)

'Identify and write Pd+1 atoms next to Al'
for item in N_list[Al]:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy.append(Atom('Pd',(zeolite[item].position[0], zeolite[item].position[1], zeolite[item].position[2]+1)))
	index = print_structure(zeolite_copy, index)

'Identify and write other oxidation states of Pd'
#to be done ...
'''

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

'''
Writing structures [one Al, two Al [NN and NNN]
'''

'single Al'
index = print_structure(zeolite, index)

'2 Al [NN]'
for item in neighbors['Si']['NN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index = print_structure(zeolite_copy, index)

'2 Al [NNN]'
for item in neighbors['Si']['NNN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index = print_structure(zeolite_copy, index)

'''
To do ...
1. verify no two structures are the same
2. verify neighbouring list always has atoms with 4 or 2 neighbors
   [check every tested Si has 4 neighbors and every tested O has 2 neighbors]
   [when I add a periodic image, am I repeated myself or identifying new combinations?]
3. Add Pd/H to NN and NNN + consider oxidation states
4. Dict to store data information

Questions
1. How can I verify two structures are not symmetric?
2. Placment of H/Pd atoms should be adjusted. How to optimize this?
3. Do metals repalce Al?
4. Should I add a criteria for distance?
5. Andrew Gettson [in Ford] did work on Al distribution on chabasize
'''
