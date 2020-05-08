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

def compare_distance(atoms1, atoms2, cutoff=0.05):
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

#def rotate_atoms(atoms, Al2, candidate):
def rotate_atoms(a, candidate):
	'''
	Generates different rotations of the structure
	Inputs:
		atoms     - ase atoms object
		Al_pos    - index of Al atom where the cell should be shifted
	Output:
		rotations - a list of atoms objected rotated around Al_index
	'''
	angle     = [0,90,180,270]
	rotations = []
	for x in angle:
		for y in angle:
			for z in angle:
				atoms = deepcopy(a)
				atoms.rotate(x,'x')
				atoms.rotate(y,'y')
				atoms.rotate(z,'z')
				rotations.append(atoms)
				#atoms.write('/home/aljama/BEA/T9/generate_structures/traj/'+candidate+'-'+str(x)+'-'+str(y)+'-'+str(z)+'.traj')

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

	output      = False 		#assume two structures do not match
	atoms1_Al, atoms2_Al = [],[] 	#index of Al atoms in the two structures

	for atom in atoms1:
		if atom.symbol == 'Al':
			atoms1_Al.append(atom.index)

	for atom in atoms2:
		if atom.symbol == 'Al':
			atoms2_Al.append(atom.index)

	for Al1 in atoms1_Al:

		Al1_xyz = deepcopy(atoms1[Al1].position)
		a1 = deepcopy(atoms1)
		a1 = shift(a1, Al1_xyz)
		#a1.write('/home/aljama/BEA/T9/generate_structures/traj/'+item)

		a2 = deepcopy(atoms2)
		rotations = deepcopy(rotate_atoms(a2, candidate))
		k = 0
		for rot in rotations:

			rot_copy = deepcopy(rot)
			for Al2 in atoms2_Al:
				rot_copy = shift(rot_copy, deepcopy(rot_copy[Al2].position))
				if candidate == '135.traj' and item == '2.traj':
					k+=1
					#rot_copy.write('/home/aljama/BEA/T9/generate_structures/traj/'+candidate+'-'+str(k)+'.traj')

				if compare_distance(a1, rot_copy) == True:
					print('overlap two strucutres!', candidate, item)
					output = True
					return output
	'''
		for Al2 in atoms2_Al:

			a2 = deepcopy(atoms2)
			rotations = deepcopy(rotate_atoms(a2, Al2, candidate))

			for rot in rotations:
				if compare_distance(a1, rot) == True:
					print('overlap two strucutres!', candidate, item)
					output = True
					return output
	'''
	return output
