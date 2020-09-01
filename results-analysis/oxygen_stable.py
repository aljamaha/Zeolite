from ase import io
from molmod import *
import os

def cation_n_O(cation_name, atoms, cutoff=2.5):
	'''
	Returns number of oxygens that are within cutoff of cation
	Inputs:
		cation_name - string of the name of the cation (assuming only one exist)
		traj_name   - name of the xyz file or traj file
		O           - list of oxygen indicies
		cutoff      - cutoff ditance between oxygen and cation
	'''

	O, n = [], 0
	atoms.write('tmp.xyz')

	'neighbor list'
	mol = Molecule.from_file('tmp.xyz')
	os.system('rm tmp.xyz')
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

	'find distance between cation and oxygen atoms neighboring Al'
	for oxygen in O:
		if atoms.get_distance(oxygen, cation_index) < cutoff:
			n+=1

	return n

