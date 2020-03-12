#!/home/mgcf/software-ws/anaconda/anaconda3/envs/molmod/bin/python

from ase import io, Atom
import os, pickle, json
from copy import deepcopy
from molmod import *
from functions import *
from qm_region import qm_region
from adsorbate import *
from check_duplicates import *

'''
Generates unique zeolite structure with 1 or 2 Al substituting Si and enumerate adsorption site
Includes options for H adsorption, Pd+1, NH3, and Pd+2
'''

'Inputs'
zeolite = io.read('../original_structures/CHA-T696.xyz')	#Zeolite structure
Pd1 	= False  #generate structures of Pd+1 (if True)
H_Z	= False	#generate structures of zeolites with H (if True)
Pd2	= False #generates sturcture of zeolites with Pd+2 (if True)
NH3	= False #genertes structures of zeolites with NH3 (if True)
Al	= 0	#index of Si atom to be replaced by an Al atom			
H_atoms = 192	#number of H atoms in original structure to account for terminal O
dir_name= 'Zeolite'

'Inputs (dont change)'
calculations    = '/home/aljama/'+dir_name+'/calculations/' #folder containing calculations
cwd  	= os.getcwd()
if os.path.exists(cwd+'/../structures_saved') == False:
	os.system('mkdir '+cwd+'/../structures_saved')
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

for item in data:
	'identify qm atoms'
	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''Detecting duplicate structures'''
data = remove_duplicates(data, struc_dir)

'''Writing structures of H-zeolites'''
if H_Z == True:
	zeolite_bare = list(data.keys())	#list of zeolites with Al but no H
	index, data  = H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list, H_atoms)

	'''identify qm region [repeated here because H in previous regions is needed for NO ads site''' 
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''Pd+2'''
if Pd2 == True:
	index, data = Pd_two(data, calculations, struc_dir, index, total_original_atoms, N_list, H_atoms)	
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''Pd+1'''
if Pd1 == True:
	index, data = Pd_one(data, struc_dir, N_list, H_atoms, index ,  total_original_atoms )
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''NH3'''
if NH3 == True:
	'''
	this needs to be adjusted
structures_so_far = list(data)

for structure in structures_so_far:

	if data[structure]['oxidation'] == 0:

		H_num = []	 #number of H atoms

		try:
			atoms = io.read(calculations+'/'+structure[0:-5]+'-opt-omegab97x-d-def2-svp/full-atoms.xyz')
		except:
			atoms    = io.read(struc_dir+'/'+structure)
			print('optimized structure for {} is not found'.format(structure))

		atoms_qm = data[structure]['qm_region']

		for atom_num in atoms_qm:
			if atoms[atom_num].symbol == 'H':
				H_num.append(atom_num)

		for H_atom in H_num:
			'prints two structures based on nearest pore in each (8 MR and 6 MR)'
			xyz_H        = atoms[H_atom].position
			H_pos        = CHA_ads(xyz_H)
			zeolite_copy = add_NH3(atoms, H_pos[2])
			index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure, adsorbate='NH3')
	'''

'''save data'''
with open(data_dir+"/data.json", "w") as write_file:
    json.dump(data, write_file, indent=4)
