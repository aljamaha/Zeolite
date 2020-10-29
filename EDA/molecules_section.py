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
mol.set_default_masses()		#reads the mass of each atom in xyz input
assert(mol.graph is None)
mol.set_default_graph()			#derive a molecular graph based on geometry

for i in mol.graph.neighbors:
		x = str(round(mol.coordinates[i][0]/angstrom,3))
		y = str(round(mol.coordinates[i][1]/angstrom,3))
		z = str(round(mol.coordinates[i][2]/angstrom,3))

		print(str(mol.symbols[i])+'\t'+x+'\t'+y+'\t'+z)

