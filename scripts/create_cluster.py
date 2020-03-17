from ase import io, Atom
from copy import deepcopy
import numpy as np
from molmod import *
import os

'''
Creates a zeolite cluster. the final structure requires using Iqmol to add missing H in Si atoms
'''

'Inputs:'
zeolite_xyz = 'original.xyz'	#name of original zeolite unit cell
x,y,z	    = 3,3,3		#copy unit cell into 3 copies in x,y,z directions
center	    = [1,1,1]		#center the cluster in this box
Trim	    = False		#further trimming of the unit cell
trim_x	    = [0.5,2.5]		#only if Trim = True
trim_y	    = [0.5,2.5]		#only if Trim = True
trim_z	    = [0.5,2.5]		#only if Trim = True

'import original trajectory file'
original_atoms  = io.read(zeolite_xyz)
atoms 		= deepcopy(original_atoms)

'get dimensions of original unit cell'
cell =  np.array(atoms.cell)

def neighbor_list(xyz_file):
	'generate a global neighbr list from atoms of the xyz_file'
	mol = Molecule.from_file(xyz_file)
	mol.set_default_masses()
	mol.set_default_graph()
	return mol.graph.neighbors

def copy_atoms(original_atoms , x , y , z):
	'''
	copy atoms to a new unit cell
	x,y,z are teh number of copies in x,y,z directions
	'''
	new_xyz = np.zeros((3))
	for index, atom in enumerate(original_atoms):
		new_xyz[0] = original_atoms[index].position[0] + cell[0][0] * x
		new_xyz[1] = original_atoms[index].position[1] + cell[1][1] * y
		new_xyz[2] = original_atoms[index].position[2] + cell[2][2] * z
		atoms.append(Atom(atom.symbol, new_xyz))

	return atoms

def clean_xyz(file_name):
	'deletes the last row in xyz file (otherwise, molmod cannot read it)'
	f = open(file_name)
	lines = f.read().split("\n")
	f.close()
	g = open('tmp.xyz','w')
	for index, line in enumerate(lines):
		if index == 0:
			g.write(line+'\n')
		elif index == 1:
			g.write('\n')
		elif index == len(lines)-1:
			continue
		else:
			g.write(line[0:-3]+'\n')
	g.close()

'shift unit cell to the center'
atoms = copy_atoms(original_atoms, center[0],center[1],center[2])

'del original atoms (since they have been shifted)'
l = []
for i in range(0,len(original_atoms)):
	l.append(i)
l.reverse()
for i in l:
	del atoms[i]

'duplicate unit cell'
for i in range(0,x):
	for ii in range(0,y):
		for iii in range(0,z):
			if i == center[0] and ii == center[1] and iii == center[2]:
				continue
			else:		
				atoms = copy_atoms(original_atoms, i,ii,iii)
atoms.write('new.xyz')


'Further trimming of the cluster'
to_be_deleted = []
if Trim == True:
	x = [cell[0][0]*trim_x[0], cell[0][0]*trim_x[1]]
	y = [cell[1][1]*trim_y[0], cell[1][1]*trim_y[1]]
	z = [cell[2][2]*trim_z[0], cell[2][2]*trim_z[1]]

	for atom in atoms:
		'z-axis'
		if atom.position[2] < trim_z[0]:
			to_be_deleted.append(atom.index)
		elif atom.position[2] > trim_z[1]:
			to_be_deleted.append(atom.index)
		elif atom.position[0] < trim_x[0]:
			'x-axis'	
			to_be_deleted.append(atom.index)
		elif atom.position[0] > trim_x[1]:
			to_be_deleted.append(atom.index)
		elif atom.position[1] < trim_y[0]:
			'y-axis'	
			to_be_deleted.append(atom.index)
		elif atom.position[1] > trim_y[1]:
			to_be_deleted.append(atom.index)
	to_be_deleted.reverse()
	for i in to_be_deleted:
		del atoms[i]
	atoms.write('trimmed.xyz')

'clean xyz so it could be read by molmod'
if Trim == True:
	clean_xyz('trimmed.xyz')
else:
	clean_xyz('new.xyz')

N_list = neighbor_list('tmp.xyz')

'Trim dangling O atoms'
O = []
for index, atom in enumerate(atoms):
	if atom.symbol == 'O':
		if len(N_list[index]) != 2:
			O.append(index)
O.reverse()
for index in O:
	del atoms[index]

'unit cell box is shifted to the center'
#atoms.cell[0][0] = cell[0][0] * (1+center[0])
#atoms.cell[1][1] = cell[1][1] * (1+center[1])
#atoms.cell[2][2] = cell[2][2] * (1+center[2])
#atoms.center(about=(center[0],center[1],center[2]))

atoms.write('missing_H.xyz')
