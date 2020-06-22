from ase import io, Atoms
import json, os

'''
prints structure of atoms nearby qm region
inputs : data dictionary + qm_structure + cutoff distance
outputs: ase atoms object
'''

cutoff = 17

'Load infomration from data dictionary (local)'
try:
	with open("dir_data.json", "r") as read_file:
	    data = json.load(read_file)
	qm = data['qm_region']
except:
	print('dir_data.json is not available')
	exit()

'load qm structure'
try:
	atoms = io.read('input.xyz')
except:
	print('input.xyz is not available')
	exit()

'generate list of atoms outside cutoff'
to_be_deleted = []
for atom in atoms:
	for qm_atom in qm:
		d = atoms.get_distance(atom.index, qm_atom)
		if d > cutoff:
			if atom.index not in to_be_deleted:
				if atom.index not in data['qm_region']:
					to_be_deleted.append(atom.index)

to_be_deleted.reverse() #so we don't mess up the numbering

def surrounding(index):

	'read input traj'
	try:
		atoms = io.read('traj-full-atoms/'+index+'.xyz')
	except:
		return 'terminate'

	'delete items'
	for item in to_be_deleted:
		del atoms[item]

	atoms.write('traj-surroundings/'+index+'.traj')

for j in range(1,1000):
	if j == 1 and os.path.exists('traj-surroundings') == False:
		os.system('mkdir traj-surroundings')
	tmp = surrounding(str(j))
	if tmp == 'terminate':
		break
	
