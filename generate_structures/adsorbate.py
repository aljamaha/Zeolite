from functions import *
from ase import io, Atom
from copy import deepcopy

def position_Pd1(Al_index, atoms, O1, O2):
	'''
	Position of Pd1 between two O atoms and at 1.7d from Al_index
	Inputs:
		Al_index - index of Al atom
		atoms	 - ase atoms object
		O1	 - index of oyxgen atom (1)
		O2	 - index of oyxgen atom (2)
	Output:
		Pd+1 position
		
	'''
	copy_atoms = deepcopy(atoms)
	mid_pt = middle_of_2_Al(atoms[O1].position, atoms[O2].position)	#mid points between two oxygen atoms
	d = mid_pt - atoms[Al_index].position	#place Pd+1 further away from Al_index
	Pd1_pos = mid_pt + 1.7*d

	return Pd1_pos


def add_metal(zeolite, ads, position, H=1.5):
	'''
	Inputs:
		zeolite   : structure of zeolite
		ads       : metal structure to be added (list format)
		position: xyz coordinates
	Outputs:
		structure of zeolite with metal adsorbed
	'''
	zeolite_copy = deepcopy(zeolite)

	if ads == 'Pd':
		zeolite_copy.append(Atom('Pd',(position[0], position[1], position[2]+H)))

	return zeolite_copy

def add_Pd1(Al_index, N_list, atoms):
	'''
	Adds Pd+1 between O1 and O2 at a 1.5d distance (where d is distance between Al and the midpoint of O1, O2)
	Inputs:
		Al_index - index of Al atom
		atoms	 - ase atoms object
		N_list   - global neighbor list of all atoms
	Outputs:
		4 copies of Structure of zeolite with Pd+1 at different positions
		
	'''
	O_N = individual_NL(Al_index, N_list)['O']['N'] #O atoms neighboring Al
	copies = []

	'Distances from O_N[0]'
	O_d = {} 					#distances from O_N[0]
	for O in O_N:
		O_d[O] = {}
		O_d[O] = atoms.get_distance(O,O_N[0])

	'Atom furthest away from O_N[0]'
	far_atom , d = '', 0
	for O in O_d:
		if O_d[O] > d:
			far_atom = O
			d = O_d[O]

	for O in O_N:
		if O != O_N[0] and O != far_atom:
			'Add Pd+1 between O_N[0] and the two Os next to it'
			Pd_pos = position_Pd1(Al_index, atoms, O_N[0], O)
			copy_atoms = deepcopy(atoms)
			copy_atoms.append(Atom('Pd', Pd_pos))
			copies.append(copy_atoms)

			'Add Pd+1 between far_atom and the two atoms (not O_N[0])'
			Pd_pos = position_Pd1(Al_index, atoms, far_atom, O)	
			copy_atoms = deepcopy(atoms)
			copy_atoms.append(Atom('Pd', Pd_pos))
			copies.append(copy_atoms)

	return copies

def Pd_H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list, H_atoms):
	'''
	write structures of Pd-H-zeolites where 2 Al atoms are present
	Inputs:
		zeolite_bare: list of zeolites without H or M
		struc_dir   : directory of structures
		data	    : dictionary of the data
		neighbors   : dictionary of the neighboring Si/Al
		index       : index of previous structure
	'''
	for item in zeolite_bare:
		atoms = io.read(struc_dir+'/'+data[item]['reference'])
		O_index = O_neighbor_indicies(atoms, N_list)
		for i in O_index[0]:
			zeolite_copy = add_H(atoms, i, 1)
			index, data = print_structure(zeolite_copy, index, data[item]['N'], item,struc_dir, data, H_atoms)
		
		for i in O_index[1]:
			zeolite_copy = add_H(atoms, i, 1)
			index, data = print_structure(zeolite_copy, index, data[item]['N'], item,struc_dir, data, H_atoms)

	return index, data

def Pd_two(data, calculations, struc_dir, index, total_original_atoms, N_list, H_atoms):
	'''
	Generates Strcutreus for Pd+2 for Al-Al pair
	Inputs:
		data         - json data file
		calculations - dir for saving calculations
		struc_dir    - dir to save structures
		index        - calculation index number
		total_original_atoms - total number of atoms in the original zeolite
		N_list       - global neighbor list
		H_atoms      - total number of H atoms in original zeolite
	Output: Generates and print structures (saves them also in data)
	'''
	structures_so_far = list(data)
	for structure in structures_so_far:

		if data[structure]['oxidation'] == -2:
		
			n_Al = 0 #number of Al atoms, position of Al atoms	
			Al_num = []	

			try:
				atoms = io.read(calculations+'/'+structure[0:-5]+'-opt-omegab97x-d-def2-svp/full-atoms.xyz')
			except:
				atoms    = io.read(struc_dir+'/'+structure)
				#print('optimized structure for {} is not found'.format(structure))

			atoms_qm = data[structure]['qm_region']	

			for atom_num in atoms_qm:
				if atoms[atom_num].symbol == 'Al':
					n_Al += 1
					Al_num.append(atom_num)

			if n_Al == 2:	
				'add to both 6 and 8 MR'
				add_Pd_pos = middle_of_2_Al(atoms[Al_num][0].position, atoms[Al_num][1].position)
				Pd_pos = CHA_ads(add_Pd_pos)
				zeolite_copy = add_metal(atoms, 'Pd', Pd_pos[0], 0) 
				index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure, metal = 'Pd')			
				zeolite_copy = add_metal(atoms, 'Pd', Pd_pos[1], 0) 
				index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure, metal = 'Pd')

	'identify qm region' 
	#for item in data:
	#	#print('='*5,'\n', item)
	#	#data = qm_region(data, item, struc_dir, N_list, total_original_atoms)

	return index, data

def H_zeolite(zeolite_bare, struc_dir, data, neighbors, index, N_list, H_atoms, Al):
	'''
	write structures of H-zeolites
	Inputs:
		zeolite_bare: list of zeolites without H or M
		struc_dir   : directory of structures
		data	    : dictionary of the data
		neighbors   : dictionary of the neighboring Si/Al
		index       : index of previous structure
	'''
	for item in zeolite_bare:
		atoms = io.read(struc_dir+'/'+str(data[item]['reference']))
		if data[item]['Al'] == 1:
			for O in neighbors['O']['N']:
				zeolite_copy = add_H(atoms, O, 1)
				index, data = print_structure(zeolite_copy, Al, index, 'N', item,struc_dir, data, H_atoms)
		else:
			O_index = O_neighbor_indicies(atoms, N_list)
			for i in O_index[0]:
				for j in O_index[1]:
					zeolite_copy = add_H(atoms, i, 2, j)
					index, data = print_structure(zeolite_copy, Al, index, data[item]['N'], item,struc_dir, data, H_atoms)

	return index, data

def add_NH3(zeolite, position):
	'''
	Inputs:
		zeolite   : structure of zeolite
		ads       : metal structure to be added (list format)
		position: xyz coordinates (preferablly returned from CHA_ads where it identified nearest pore)
	Outputs:
		structure of zeolite with NH3 adsorbed
	'''
	zeolite_copy = deepcopy(zeolite)
	NH = 1	#distance between N-H

	zeolite_copy.append(Atom('N',(position[0]    , position[1]    , position[2])))
	zeolite_copy.append(Atom('H',(position[0]+NH, position[1]    , position[2])))
	zeolite_copy.append(Atom('H',(position[0]-NH, position[1]    , position[2])))
	zeolite_copy.append(Atom('H',(position[0]    , position[1]+NH, position[2])))

	return zeolite_copy

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

def Pd_one(data, struc_dir, N_list, H_atoms, index , total_original_atoms ):
	structures_so_far = list(data)
	for structure in structures_so_far:
		if data[structure]['oxidation'] == -1:	
			'add only Pd+1'
			atoms    = io.read(struc_dir+'/'+structure)
			Al_index = []
			for atom in atoms:
				'find Al atom'
				if atom.symbol == 'Al':
					Al_index = atom.index
					break
			zeolite_copies = add_Pd1(Al_index, N_list, atoms)
			for zeolite_copy in zeolite_copies:
				index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure)
				
		elif  data[structure]['oxidation'] == -2:
			'add Pd+1 and H'
			atoms    = io.read(struc_dir+'/'+structure)
			Al_index = []
			for atom in atoms:
				'find Al atom'
				if atom.symbol == 'Al':
					Al_index.append(atom.index)

			for j, Al_atom in enumerate(Al_index):
				'loop over each of the two Al atoms'
				for O_n in individual_NL(Al_atom, N_list)['O']['N']:
					'loop over each O neighboring the Al atom, add H to it'
					zeolite_copy = add_H(atoms, O_n, 1)	
					'Add Pd on the opposite Al atom'
					if j == 0:
						zeolite_copies = add_Pd1(Al_index[1], N_list, zeolite_copy)
					else:
						zeolite_copies = add_Pd1(Al_index[0], N_list, zeolite_copy)

					for zeolite_copy in zeolite_copies:
						index, data  = print_structure(zeolite_copy, index, data[structure]['N'], data[structure]['reference'] , struc_dir, data, H_atoms, reference_H = structure)

	'identify qm region [repeated here because H in previous regions is needed for NO ads site'
	#for item in data:
	#	data = qm_region(data, item, struc_dir, N_list, total_original_atoms)
		
	return index, data


def NH3_old():
	'''
	'the original NH3 function in main script'
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
