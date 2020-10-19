#import packages
import sys, json
from molmod import *	#further info are here: http://molmod.github.io/molmod/tutorial/install.html
from ase import io

'''
Objective: prepare external charges section in EDA calc. for QM/MM calculations
'''

'Inputs'
input_xyz = 'qm-final.xyz'

'qm reigon (Assumes dir_data.json file is available)'
mol = Molecule.from_file(input_xyz)	#reads the input <xyz coordinates>
#with open('dir_data.json','r') as read:
#	data = json.load(read)
#	qm_atoms  = data['qm_region']

mol.set_default_masses()		#reads the mass of each atom in xyz input
assert(mol.graph is None)
mol.set_default_graph()			#derive a molecular graph based on geometry

'''
'Identiy O-atoms (in MM region) connected to qm atoms [replaced by H for the QM calculation]'
link_O = []
for i in mol.graph.neighbors:
	'loop over all atoms in xyz input file'
	if i in qm_atoms:
		indexes = [n for n in mol.graph.neighbors[i]]
		for j in range(0,len(indexes)):
        	        if not indexes[j] in qm_atoms:
                	   link_O.append(indexes[j])

'Identify Si atoms attached to link_O in MM region'
'[Si atoms connected to linking O atoms (atom type -23/-22/-21 depending on number of linking O atoms the Si atom is connected to)]'
link_Si = []
for i in mol.graph.neighbors:
	'loop over all atoms in xyz input file'
	if i in link_O:
		indexes = [n for n in mol.graph.neighbors[i]]
		for j in range(0,len(indexes)):
			if not indexes[j] in qm_atoms and mol.symbols[indexes[j]] == 'Si':
				link_Si.append(indexes[j])
'''
for i in mol.graph.neighbors:
		#indexes = [n for n in mol.graph.neighbors[i]]
		#if i in qm_atoms:
		#type = 0 #no MM information is needed here

		x = str(round(mol.coordinates[i][0]/angstrom,3))
		y = str(round(mol.coordinates[i][1]/angstrom,3))
		z = str(round(mol.coordinates[i][2]/angstrom,3))

		print(str(mol.symbols[i])+'\t'+x+'\t'+y+'\t'+z)

