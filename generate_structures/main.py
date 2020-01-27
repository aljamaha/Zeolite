#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io, Atom
import os, pickle, json
from copy import deepcopy
from molmod import *
from functions import *
from qm_region import qm_region

'''
Generates unique zeolite structure with 1 or 2 Al substituting Si and enumerate adsorption sites [H and metal]
'''

'Inputs'
zeolite = io.read('../original_structures/CHA-T696.xyz')	#Zeolite structure
Al	= 0	#index of Si atom to be replaced by an Al atom			
H_atoms = 192	#number of H atoms in original structure to account for terminal O
metals = {}
metals['PdO']   = {'composition':['Pd','O'], 'oxidation_state':[-2,0,2]}
metals['Pd2O']  = {'composition':['Pd','O','Pd'], 'oxidation_state':[-2,2,6]}
metals['PdO2']  = {'composition':['Pd','O','O'], 'oxidation_state':[-4,-2,0]}
metals['Pd2O2'] = {'composition':['Pd','O','Pd','O'], 'oxidation_state':[-4,0,4]}
metals['Pd']    = {'composition':['Pd'], 'oxidation_state':[0,2,4]}
calculations    = '/home/aljama/CHA/calculations/' #folder containing calculations

'Inputs (dont change)'
cwd  	= os.getcwd()
struc_dir = cwd+'/../structures_saved'	#dir to store structures
index 	= 0			#index of the structure
data 	= {}			#store details of each structure
neighbors = {}			#storing neighbors for Si and O
neighbors['O']  = {'N':[],'NN':[],'NNN':[]}
neighbors['Si'] = {'N':[],'NN':[],'NNN':[]}
data_dir        = cwd+'/../data'
total_original_atoms = len(zeolite) 	#number of total atoms in the zeolite structure

'substitute Si with Al'
zeolite[Al].symbol = 'Al'

'neighbour list'
zeolite.write('tmp.xyz')
N_list = neighbor_list('tmp.xyz')	#dict of neighbors list
os.system('rm tmp.xyz')

'''Building Si and O [N, NN, and NNN]'''
neighbors = identify_N(N_list[Al],N_list, neighbors, Al)
neighbors = identify_NN_O(neighbors, N_list)
neighbors = identify_NN_Si(neighbors, N_list)
neighbors = identify_NNN_O(N_list, neighbors)
neighbors = identify_NNN_Si(N_list, neighbors)

'''Writing structures [one Al, two Al][NN and NNN]'''
'single Al'
index, data = print_structure(zeolite, index, 'N', str(index+1)+'.traj' , struc_dir, data, H_atoms)

'2 Al [NN]'
for item in neighbors['Si']['NN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index, data = print_structure(zeolite_copy, index, 'NN', str(index+1)+'.traj',struc_dir, data, H_atoms)

'2 Al [NNN]'
for item in neighbors['Si']['NNN']:
	zeolite_copy = deepcopy(zeolite)
	zeolite_copy[item].symbol = 'Al'
	index, data = print_structure(zeolite_copy, index, 'NNN', str(index+1)+'.traj',struc_dir, data, H_atoms)

'''Writing structures of H-zeolites'''
zeolite_bare = list(data.keys())	#list of zeolites with Al but no H
index, data  = H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list, H_atoms)

'''identify qm region [repeated here because H in previous regions is needed for NO ads site''' 
for item in data:
	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

''''writing structures of NH3 adsorbed on H'''
structures_so_far = list(data)
for structure in structures_so_far:
	if data[structure]['oxidation'] == 0:
		try:
			atoms = io.read(calculations+'/'+structure[0:-5]+'-opt-omegab97x-d-def2-svp/full-atoms.xyz')
		except:
			atoms    = io.read(struc_dir+'/'+structure)
		atoms_qm = data[structure]['qm_region']	
		for atom_num in atoms_qm:
			atom_type = atoms[atom_num].symbol
			if atom_type == 'H':	
	
				zeolite_copy = add_ads(atoms, ['N','H','H','H'], atoms[atom_num].position, index) 
				index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure, adsorbate='NH3')	

'''identify qm region [repeated here because H in previous regions is needed for NO ads site''' 
for item in data:
	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''save data'''
with open(data_dir+"/data.json", "w") as write_file:
    json.dump(data, write_file, indent=4)

print('exited after NH3 ads ..')

exit()

'''Writing structures of metal modified zeolites'''
no_metal_zeolite = list(data) #List of structures with no introduced metal [includes ones with H]

for structure in no_metal_zeolite:
	if data[structure]['oxidation'] == -2:
		'oxidation stae of +2 is needed'
		'do I need to do other oxidation states?'
		for comp in metals:
			for ox in metals[comp]['oxidation_state']:
				if ox == 2:
					atoms = io.read(struc_dir+'/'+structure)
					for atom in atoms:
						if atom.symbol == 'Al':	
							zeolite_copy = add_metal(atoms, metals[comp]['composition'], atom.position)
							index, data = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'], struc_dir, data, H_atoms)

'''identify qm region'''
for item in data:
	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''NO adsorption'''
no_ads_zeolites = list(data)
for structure in no_ads_zeolites:
	if data[structure]['oxidation'] == 0:
		atoms = io.read(struc_dir+'/'+structure)
		atoms_qm = data[structure]['qm_region']
		for atom_num in atoms_qm:
			atom_type = atoms[atom_num].symbol
			if atom_type == 'H':
				zeolite_copy = add_ads(atoms, ['N','O'], atoms[atom_num].position)
				index, data = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure)

'''identify qm region'''
for item in data:
	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''save data'''
with open(data_dir+"/data.json", "w") as write_file:
    json.dump(data, write_file, indent=4)

'''
Questions:
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
