from ase import io, Atoms
import json

'Reorder atoms that QM atoms at the beginning and MM atoms at the end'

'create atoms object to save new data'
ordered_atoms  = Atoms('Au', positions=[[0, 0, 0]])

with open("dir_data.json", "r") as read_file:
    data = json.load(read_file)

atoms = io.read('input.xyz')

'new atoms numbers'
d = 0
new_qm = []
for i in range(0,len(atoms)):
	if i in data['qm_region']:
		ordered_atoms.append(atoms[i])
		new_qm.append(d)	
		d += 1

for i in range(0,len(atoms)):
	if i not in data['qm_region']:
		ordered_atoms.append(atoms[i])
		d += 1

data['qm_region'] = new_qm

with open("dir_data_qmmm_order.json", "w") as write_file:
    json.dump(data, write_file, indent=4)

'del Au atom'
del ordered_atoms[0]

'write ordered qm atoms'
ordered_atoms.write('input-ordered.xyz')
