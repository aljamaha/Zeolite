from ase import io, Atoms
import json, sys

'''
Generates structures centered on T atom (to visualize local environment)
'''

cutoff = 7 #distance greater then this, all atoms will be excluded
T_atom = [128,136,144,152,160,168,176,184,188, 163,149,191,153,170]

for T in T_atom:
	atoms = io.read('3x3x3.xyz')
	#atoms = io.read('BEA.cif')
	atoms[T].symbol = 'Pd'

	'generate list of atoms outside cutoff'
	to_be_deleted = []
	for atom in atoms:
		d = atoms.get_distance(atom.index,T)
		if d > cutoff:
			if atom.index not in to_be_deleted:
				to_be_deleted.append(atom.index)			
	to_be_deleted.reverse() #so we don't mess up the numbering

	'delete items'
	for item in to_be_deleted:
		del atoms[item]

	atoms.write('T-'+str(T)+'.traj')
