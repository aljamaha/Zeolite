from ase import io, Atoms
import json
from molmod import *
from functions import *

'Description: generate structures of few added qm atoms based on a previous structure'

'Inputs'
with open('dir_data.json','r') as read:
	data = json.load(read)
qm_region = data['qm_region']

'General Inputs'
NL    = neighbor_list('input.xyz')
atoms = io.read('input.xyz')
candidates = []

'Find potential Si atoms connected to qm region'
for atom in atoms:
	if atom.index in qm_region:
		tmp, candidates = False, []
		if atom.symbol in ['Al','Si']:
			Si_N = individual_NL(atom.index, NL)['Si']['N']
			for index in Si_N:
				if index not in qm_region:
					tmp = True
					candidates.append(index)
			if tmp == True:
				#if len(candidates) <3:
				print(atom.index, atom.symbol, atom.position, candidates)

'identify O atoms connected to Si in the qm region'
for item in NL:
		if set(NL[item]) <= set(qm_region):
			if item in qm_region:
				continue
			else:
				qm_region.append(item)
'write new unit cell'
new   = Atoms('Ni', [(0, 0, 0)],cell=[1, 1, 1]) #dummy structure
for j in qm_region:
	new.append(atoms[j])
del new[0]

new.write('qm-initial.traj')
print('Length of qm region:', len(new))

data['qm_region'] = qm_region
with open('dir_data.json', 'w') as write:
	json.dump(data, write, indent=4)
