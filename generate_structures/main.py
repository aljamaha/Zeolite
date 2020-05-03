#!/home/mgcf/software-ws/anaconda/anaconda3/envs/molmod/bin/python

from ase import io, Atom
import os, pickle, json
from copy import deepcopy
from molmod import *
from functions import *
from qm_region import qm_region
from adsorbate import *
from check_duplicates import *
import time

'''
Generates unique zeolite structure with 1 or 2 Al substituting Si and enumerate adsorption site
Includes options for H adsorption, Pd+1, NH3, and Pd+2
'''
start = time.time()
'Inputs'
zeolite_original = io.read('../original_structures/original.xyz')	#Zeolite structure
Pd1 	= False  #generate structures of Pd+1 (if True)
H_Z	= True	#generate structures of zeolites with H (if True)
Pd2	= False #generates sturcture of zeolites with Pd+2 (if True)
T_atom	= [128,136,144,152,160,168,176,184,188]	#index of Si atom to be replaced by an Al atom	
T_atom	= [128]	#index of Si atom to be replaced by an Al atom	
H_atoms = 348	#number of H atoms in original structure to account for terminal O
dir_name= 'BEA/H'
cutoff  = 0	#adds O atoms < cutoff distance to qm region
n_MR_max = 7	#maximum number of MR of interest
n_atoms_original = 192	#number of atoms in a unit cell

'Inputs (dont change)'
calculations    = '/home/aljama/'+dir_name+'/calculations/' #folder containing calculations
cwd  	= os.getcwd()
if os.path.exists(cwd+'/../structures_saved') == False:
	os.system('mkdir '+cwd+'/../structures_saved')
struc_dir = cwd+'/../structures_saved'	#dir to store structures
index 	= 0			#index of the structure
data 	= {}			#store details of each structure
data_dir        = cwd+'/../data'
total_original_atoms = len(zeolite_original) 	#number of total atoms in the zeolite structure

for Al in T_atom:

	print('T-atom:', Al)
	'initialize neighbors dict'
	
	zeolite = deepcopy(zeolite_original)
	
	'substitute Si with Al'
	zeolite[Al].symbol = 'Al'

	'neighbour list'
	zeolite.write('tmp.xyz')
	N_list = neighbor_list('tmp.xyz')	#dict of neighbors list
	os.system('rm tmp.xyz')

	'Building Si and O [N, NN, and NNN]'
	neighbors   = individual_NL(Al, N_list)

	'Writing structures [one Al, two Al][NN and NNN]'
	'single Al'
	index, data = print_structure(zeolite, Al, index, 'N', str(index+1)+'.traj' , struc_dir, data, H_atoms, Al)

	'2 Al [NN]'
	for item in neighbors['Si']['NN']:
		zeolite_copy = deepcopy(zeolite)
		zeolite_copy[item].symbol = 'Al'
		index, data = print_structure(zeolite_copy, Al, index, 'NN', str(index+1)+'.traj',struc_dir, data, H_atoms, Al)

	'2 Al [NNN]'
	for item in neighbors['Si']['NNN']:
		zeolite_copy = deepcopy(zeolite)
		zeolite_copy[item].symbol = 'Al'
		index, data = print_structure(zeolite_copy, Al, index, 'NNN', str(index+1)+'.traj',struc_dir, data, H_atoms, Al)


end = time.time() - start
print('Time to generate structures ..', end)

'removes structures with Al-Al as neighbors'
data = Al_Al_N(struc_dir, data, N_list)

'Limit to only single unit cell'
#data = unit_cell_limit(data, struc_dir, n_atoms_original)

'identify qm atoms, and add Al-Al distance to data.json'
print('identifying qm region .. ')
for item in data:
	data  = qm_region(data, item, struc_dir, N_list, total_original_atoms,cutoff,n_MR_max, turn_cutoff='off')
	atoms = io.read(struc_dir+'/'+item) 
	data[item]['Al-Al distance'] = Al_Al_distance(atoms)

print('Time to identify qm region',  time.time() - start)

'''Detecting duplicate structures'''
print('Deleting duplicates ..')
data = remove_duplicates(data, struc_dir)

'''Writing structures of H-zeolites'''
if H_Z == True:
	print('creating H structures ...')
	zeolite_bare = list(data.keys())	#list of zeolites with Al but no H
	index, data  = H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list, H_atoms, Al)
	#index, data  = H_zeolite(zeolite_bare, Al, struc_dir, data, neighbors, index, N_list, H_atoms)

	'''identify qm region [repeated here because H in previous regions is needed for NO ads site''' 
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms, cutoff, n_MR_max, turn_cutoff='on')

print('Time for H atoms ...',  time.time() - start)

'''Pd+2'''
if Pd2 == True:
	print('creating Pd+2 structures ...')
	index, data = Pd_two(data, calculations, struc_dir, index, total_original_atoms, N_list, H_atoms)	
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''Pd+1'''
if Pd1 == True:
	print('creating Pd+1 structures ...')
	index, data = Pd_one(data, struc_dir, N_list, H_atoms, index ,  total_original_atoms )
	for item in data:
		data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

'''save data'''
with open(data_dir+"/data.json", "w") as write_file:
    json.dump(data, write_file, indent=4)
