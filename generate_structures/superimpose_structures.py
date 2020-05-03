from ase import io
import os
from copy import deepcopy
from ase import Atoms

def shift(atoms, Al_position):
	'''
	shifts the structure to the origin (where origin is Al_position)
	Inputs:
		atoms       - ase atoms object
		Al_position - position of Al atom where the unit cell will be centered
	'''
	for i in atoms:
		'moves each atom so the origin is Al_position'
		i.position = i.position - Al_position
	return atoms

def compare_distance(atoms1, atoms2, cutoff=0.01):
	'''
	checks if two structures are identical by comparing distances of each atom
	Inputs:
		atoms1/atoms2 - ase object for each atom to be compared
		cutoff        - threshold for comparing the distances between thet wo atoms
	Output:
		True  -  two structures are identical
		False -  two structures are not identical
	'''
	tmp   = True 	#if True, they are identical
	for atom1 in atoms1:
		d_tmp = False   #does it pass the cutoff distance test
		for atom2 in atoms2:
			d = abs(atom1.position - atom2.position)
			if d[0] < cutoff and d[1] < cutoff and d[2] < cutoff:
				'are the two atoms in xyz positions identical?'
				d_tmp = True
		if d_tmp == False:
			'verify that for each atom1, there is an atom in atoms2 that is at identical position'
			tmp = False
			break
	return tmp

def rotate(atoms, Al_pos):
	'''
	Generates different rotations of the structure
	Inputs:
		atoms     - ase atoms object
		Al_pos    - pos of Al atom where the cell should be shifted
	Output:
		rotations - a list of atoms objected rotated around Al_index
	'''
	angle     = [0,90,180,270]
	rotations = []
	for x in angle:
		for y in angle:
			for z in angle:
				atoms.rotate(x,'x')
				atoms.rotate(y,'y' )
				atoms.rotate(z,'z')
				atoms = shift(atoms, Al_pos)
				rotations.append(atoms)

	return rotations

def superimpose_structures(item, candidate, struc_dir, data):
	'''
	superimpose two structures on each other by generating many rotations on structure 2 (candidate)
	Inputs:
		item      - name of structure in struc_dir
		candidate - name of structure in struc_dir
	Output:
		True  - the two structures overlap
		False - the two structures do not overlap 
	'''
	
	'atoms1 qm region'
	full_atoms1  = io.read(struc_dir+'/'+item)
	atoms1       = Atoms('Ni', [(0, 0, 0)],cell=[1, 1, 1]) #dummy structure
	for qm in data[item]['qm_region']:
		atoms1.append(full_atoms1[qm])
	del atoms1[0]

	'atoms2 qm region'
	full_atoms2  = io.read(struc_dir+'/'+candidate)
	atoms2       = Atoms('Ni', [(0, 0, 0)],cell=[1, 1, 1]) #dummy structure
	for qm in data[candidate]['qm_region']:
		atoms2.append(full_atoms2[qm])
	del atoms2[0]

	os.chdir('/home/aljama/BEA/H/generate_structures/')
	atoms1.write(item+'.traj')
	atoms2.write(candidate+'.traj')
	exit()

	'''
	atoms1      = io.read(struc_dir+'/'+item)
	atoms1_copy = deepcopy(atoms1)
	atoms1_qm    = deepcopy(data[item]['qm_region'])
	atoms1_qm.sort(reverse = True) 
	print(atoms1_copy)
	for atom in atoms1_copy:
		if atom.index not in atoms1_qm:
			del atoms1[atom.index]


	'atoms2 qm region'
	atoms2      = io.read(struc_dir+'/'+candidate)	
	atoms2_copy = deepcopy(atoms2)
	atoms2_qm    = deepcopy(data[candidate]['qm_region'])
	atoms2_qm.sort(reverse = True) 
	for atom in atoms2_copy:
		if atom.index not in atoms2_qm:
			del atoms2[atom.index]
	'''

	output      = False 		#assume two structures do not match
	atoms1_Al, atoms2_Al = [],[] 	#index of Al atoms in the two structures

	for atom in atoms1:
		if atom.symbol == 'Al':
			atoms1_Al.append(atom.index)

	for atom in atoms2:
		if atom.symbol == 'Al':
			atoms2_Al.append(atom.index)

	for Al1 in atoms1_Al:

		Al1_xyz = atoms1[Al1].position
		a1 = deepcopy(atoms1)
		a1 = shift(a1, Al1_xyz)

		for Al2 in atoms2_Al:

			Al2_xyz = atoms2[Al2].position
			a2 = deepcopy(atoms2)
			a2 = shift(a2, Al2_xyz)
			rotations = deepcopy(rotate(a2, Al2_xyz))

			for ind, rot in enumerate(rotations):
				if compare_distance(a1, rot) == True:
					print('overlap two strucutres!', candidate, item)
					output = True
					return output
	return output
