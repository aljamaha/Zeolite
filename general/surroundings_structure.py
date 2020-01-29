from ase import io, Atoms
import json, sys

'''
prints structure of atoms within cutoff distance of first Al atom (assuming atom 1 is Al)
inputs : 
	 cutoff distance
	 [optional] xyz corrdinates of traj file - otherwise default is input.xyz
outputs: ase atoms object
'''

cutoff = 15 #distance greater then this, all atoms will be excluded

try:
	len(sys.argv[1])
	atoms = io.read('full-atoms.xyz')
except:
	'read input traj'
	try:
		atoms = io.read('input.xyz')
	except:
		atoms = io.read('input.traj')
		atoms.write('input.xyz')

'generate list of atoms outside cutoff'
to_be_deleted = []
for atom in atoms:
	#for qm_atom in qm:
	d = atoms.get_distance(atom.index,1)
	if d > cutoff:
		if atom.index not in to_be_deleted:
			to_be_deleted.append(atom.index)			
to_be_deleted.reverse() #so we don't mess up the numbering

'delete items'
for item in to_be_deleted:
	del atoms[item]

try:
	len(sys.argv[1])
	atoms.write('full-atoms_surroundings.traj')
except:
	atoms.write('surroundings.traj')

