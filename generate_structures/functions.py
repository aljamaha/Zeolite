#!/Users/hassanaljama/opt/anaconda3/bin/python

from ase import io, Atom
import os
from copy import deepcopy
from molmod import *
import pickle
from ase.build import add_adsorbate

def individual_NL(index, N_list):
	'creates a neighboring list for atom of interest (needs only atom index and complete NL'

	n_list = {} #local neighbor list specific to Si atoms next to O next to terminal Si
	n_list['Si'], n_list['O'] = {},{}
	n_list['Si']['N'], n_list['O']['N'] = [],[]
	n_list['Si']['NN'], n_list['O']['NN'] = [],[]
	n_list['Si']['NNN'], n_list['O']['NNN'] = [],[]
	n_list = identify_N(N_list[index], N_list, n_list, index)
	n_list = identify_NN_O(n_list, N_list )
	n_list = identify_NN_Si(n_list, N_list)
	n_list = identify_NNN_O(N_list, n_list)
	n_list = identify_NNN_Si(N_list, n_list)

	return n_list

def neighbor_list(xyz_file):
	'''
	Inputs:  xyz coordinates
	Outputs: dictionary of the neighbor list (for all atoms)
		e.g. [atom index]: [neighboring atoms]
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

def print_structure(atoms, Al, index, N, reference, struc_dir, data, H_atoms,reference_H=[], adsorbate='', metal=''):
	'''
	Inputs:
		atoms: ase atoms object
		index: the index of the previous structure
		N    : N, NN, or NNN
		reference: original structure without modification
		struc_dir: directory to store structures
		data   : json data file 
		H_atoms: number of terminal H atoms in original zeolite 
		adsorbate     : name of the adsorbate
	Outputs:
		prints the traj file in struc_dir folder
		returns the index of the traj file in struc_dir and data dictionary
		appends details of the structure to data dictionary
	'''

	index += 1
	name = str(index)+'.traj'
	os.chdir(struc_dir), atoms.write(name)
	data[name] = count_elements(atoms)
	data[name]['N'] = N
	data[name]['reference'] = reference
	if adsorbate != '':
		data[name]['adsorbate'] = adsorbate
	if 'H' not in list(data[name]):
		data[name]['H'] = 0
	if 'Pd' in data[name]:
		data[name]['oxidation'] = 0
	elif adsorbate != '':
		data[name]['oxidation'] = 0
		data[name]['reference_H'] = reference_H
	else:	
		data[name]['oxidation'] = int(data[name]['Al'])*3 - H_atoms + (int(data[name]['H']) - H_atoms) + int(data[name]['Si'])*4 - int(data[name]['O'])*2 
	data[name]['total_atoms'] = len(atoms)
	if metal != '':
		data[name]['metal'] = metal
	data[name]['T-site'] = Al
		

	return index, data

def distance_to_others(atoms, input_atom_index, cutoff = 0.8):
	'''
	finds the distance of the atom wrt other atoms
	Inputs  : input_atom_index: index of atom of interest
		  atoms: ase atoms object
	Outputs : regular (no atom too close)
		  close   (some atoms might be too close)
	
	'''

	d = 10 #large initial value
	for atom in atoms:
		if atom.index != input_atom_index:
			if atom.symbol != 'H':
				d = atoms.get_distance(atom.index, input_atom_index)
				if abs(d) < cutoff:
					distance = d
					break
	return d
	

def O_neighbor_indicies(atoms, N_list):
	'''
	Inputs:  ase atoms object, neighbor list (N_list)
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
		#print(i)
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

	return round(distance,3)

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

def middle_of_2_Al(atom1, atom2):
	'''
	finds the middle distance between two Al atoms
	inputs:
		atom1: position of first Al atom
		atom2: position of second Al atom
	output: 
		position to add ads
	'''
	
	pos = (atom1+atom2)/2.0
	
	return pos
	
def CHA_ads(xyz):
	'''
	find the nearest pore for adsorbing an atom on Chabasite [6 MR and 8 MR]
	Inputs:
		xyz - coordinates of the atoms to be adsorbed
	Output:
		pore = xyz of nearest pore in 6 and 8 MR, and the shorted of the two
	'''

	#dict of compiled xyz postions where pore is located (either along y or z axis)
	CHA = {}
	CHA['y'] = [[-2,0, -5], [10,0, -5], [6,0,0], [-6,1.5, 0], [2,0,5], [10,0,10], [-2,0,-5], [-2,-2,10]]
	CHA['z'] = [[0,0,1.5], [4,7, 3], [4,-7, 3], [8,0,8], [-8,0,0], [-4,7,0], [-4,-7,3], [12,-7,3]]

	d1, d2 = 10, 10 #assign a large initial value
	pore1, pore2 = [100,100,100], [100,100,100]
	for item in CHA['z']:
		distance = np.sqrt( (item[0] - xyz[0])**2 + (item[1] - xyz[1])**2 )
		if distance < d1:
			d1 = distance
			pore1 = item

	for item in CHA['y']:
		distance = np.sqrt( (item[0] - xyz[0])**2 + (item[2] - xyz[2])**2 )
		if distance < d2:
			d2 = distance
			pore2 = item

	#retain position in the original axis
	pore1 = [pore1[0], pore1[1], xyz[2]]
	pore2 = [pore2[0], xyz[1]  , pore2[2]]

	if d2>d2:
		pore_min = pore1
	else:
		pore_min = pore2

	return pore1, pore2, pore_min
				
def O_cutoff(Al_indicies, data, traj, atoms, cutoff):
	'''
	identifies O atoms within a cutoff from each Al	
	Inputs:
		Al_indicies - indexes of Al atoms
		data	    - global json data list
		traj	    - name of traj file
		atoms	    - ase atoms object
		cutoff	    - cutoff beyond which atoms are not included
	Output:
		Updated data[traj][qm_region] with O-atoms <cutoff distance from Al atoms now included
	'''
	for Al in Al_indicies:
		for atom in atoms:
			if atom.symbol == 'O':
				d = atoms.get_distance(Al,atom.index)
				if d < cutoff:
					if atom.index not in data[traj]['qm_region']:	
						#print(atom.index)
						data[traj]['qm_region'].append(atom.index)
	return data

def Al_Al_N(struc_dir, data, N_list):
	'''
	find Al-Al pairs that are neighbors and deletes them
	Inputs:
		item - structure traj name
		data - global json data file
	Output: return json data with deleted entries of Al-Al pairs
	'''

	copy_data = deepcopy(data)	#so I could iterate over it, and t

	for item in copy_data:
		tmp = False
		if data[item]['Al'] == 2:
			atoms = io.read(struc_dir+'/'+item)
			for atom in atoms:
				if atom.symbol == 'Al':
					Al_index = atom.index
					break
			n_Al = individual_NL(Al_index, N_list)
			for index in n_Al['Si']['N']:
				if atoms[index].symbol == 'Al':
					tmp = True
		if tmp == True:
			del data[item]

	return data
	
def unit_cell_limit(data, struc_dir, n_atoms_original):
	'''
	limits substituted Al to be only in the same unit cell
	Inputs:
		data - json data files
		struc_dir - dir where structures are saved
		n_atoms_original - number of atoms in the original unit cell
	'''
	copy_data = deepcopy(data)
	for item in copy_data:
		atoms = io.read(struc_dir+'/'+item)
		for atom in atoms:
			if atom.symbol == 'Al' and atom.index > n_atoms_original:
				try:
					del data[item]
					print(item, 'not in same unit cell', atom.index)
				except:
					print('could not be deleted')
	return data
