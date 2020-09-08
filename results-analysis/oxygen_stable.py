from ase import io
from molmod import *
import os

def cation_n_O(cation_name, atoms, cutoff1=2.5, cutoff2=4.4):
	'''
	Returns number of oxygens that are within cutoff of cation 
	Inputs:
		cation_name - string of the name of the cation (assuming only one exist)
		atoms       - ase atoms object
		cutoff1      - cutoff ditance between oxygen and cation
		cutoff2      - cutoff ditance between oxygen and cation
	Output:
		O_Si	    - # of oxygen atoms next to Si (<cutoff2)
		O_Al	    - # of oxygen atoms next to Al (<cutoff1)
		distances   - distance of oxygen atoms to cation
	'''


	O, O_Al, O_Si, distances =  [], 0,0,[]
	atoms.write('tmp.xyz')

	'neighbor list'
	mol = Molecule.from_file('tmp.xyz')
	mol.set_default_masses()
	mol.set_default_graph()
	N_list = mol.graph.neighbors

	'find Al atoms and their neighboring oxygens'
	for atom in atoms:
		if atom.symbol == 'Al':
			tmp = N_list[atom.index]
			for item in tmp:
				O.append(item)

	'find cation'
	for atom in atoms:
		if atom.symbol == cation_name:
			cation_index = atom.index
			break

	'find oxygens next to Si and Al (within the first cutoff)'
	for atom in atoms:
		if atom.symbol == 'O':
			d = atoms.get_distance(atom.index, cation_index) 
			if d < cutoff1:
				distances.append(round(d,1))
				if atom.index in O:
					O_Al+=1
				else:
					O_Si+=1

	'find oxygens next to Si and Al (between cutoff1 and cutoff2'
	for atom in atoms:
		if atom.symbol == 'O':
			d = atoms.get_distance(atom.index, cation_index) 
			if d > cutoff1 and d < cutoff2:
				distances.append(round(d,1))
				O_Si+=1

	return O_Al, O_Si, distances
