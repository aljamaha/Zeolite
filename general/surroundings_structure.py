from ase import io, Atoms
import json

'''
prints structure of atoms nearby qm region
inputs : data dictionary + qm_structure + cutoff distance
outputs: ase atoms object
'''

cutoff = 15 #distance greater then this, all atoms will be excluded

'Load infomration from data dictionary (local)'
try:
	with open("dir_data.json", "r") as read_file:
	    data = json.load(read_file)
except:
	print('dir_data.json is not available')
	exit()

'load qm structure'
try:
	qm  = io.read('qm-initial.traj')
except:
	print('qm-initial.traj is not available')
	exit()

'read input traj'
try:
	atoms = io.read('input.xyz')
except:
	atoms = io.read('input.traj')
	atoms.write('input.xyz')

'generate list of atoms outside cutoff'
to_be_deleted = []
for atom in atoms:
	for qm_atom in qm:
		d = atoms.get_distance(atom.index,qm_atom.index)
		if d > cutoff:
			if atom.index not in to_be_deleted:
				to_be_deleted.append(atom.index)			
to_be_deleted.reverse() #so we don't mess up the numbering

'delete items'
for item in to_be_deleted:
	del atoms[item]

atoms.write('surroundings.traj')

