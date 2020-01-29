import os, sys
from molmod import *
from ase import io

'''
Objective:    convert qchem output to traj files that only includes atoms surrounding first Al atom
Inputs:	      cutoff of atoms surrounding first Al atom
Output:       xyz files under the folder surroundings-traj
'''

cutoff = 10 #distance greater then this, all atoms will be excluded
to_be_deleted = []

if os.path.exists('traj-surrounding-atoms') == False:
	os.system('mkdir traj-surrounding-atoms')	#data are saved here

for i in range(1,10000):
	if os.path.isfile('traj-full-atoms/'+str(i)+'.xyz') == True:
		atoms = io.read('traj-full-atoms/'+str(i)+'.xyz')
	else:
		break
		
	if to_be_deleted == []:
		'generate list of atoms outside cutoff'
		for atom in atoms:
			d = atoms.get_distance(atom.index,1)
			if d > cutoff:
				if atom.index not in to_be_deleted:
					to_be_deleted.append(atom.index)
		to_be_deleted.reverse() #so we don't mess up the numbering

	'delete items'
	for item in to_be_deleted:
		del atoms[item]

	atoms.write('traj-surrounding-atoms/'+str(i)+'.traj')

