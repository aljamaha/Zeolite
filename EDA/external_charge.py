#import packages
import sys, json
from molmod import *	#further info are here: http://molmod.github.io/molmod/tutorial/install.html
from ase import io

'''
Objective: prepare external charges section in EDA calc. for QM/MM calculations
'''

'Inputs'
input_xyz = 'input.xyz'		#full Q/MM

'qm region (assumes json file is available)'
mol = Molecule.from_file(input_xyz)	#reads the input <xyz coordinates>
with open('dir_data.json','r') as read:
	data = json.load(read)
	qm_atoms  = data['qm_region']

'Defaults'
atoms = io.read(input_xyz)
mol.set_default_masses()		#reads the mass of each atom in xyz input
assert(mol.graph is None)
mol.set_default_graph()			#derive a molecular graph based on geometry

'charge dict'
charge_dict = {}
charge_dict[-1]  = -0.35     
charge_dict[-2]  = 0.7
charge_dict[-3]  = -0.175   
charge_dict[-21] =  0.175
charge_dict[-22] =  0.35
charge_dict[-23] =  0.525

'Identiy O-atoms (in MM region) connected to qm atoms'
link_O = []
for i in mol.graph.neighbors:
	'loop over all atoms in xyz input file'
	if i in qm_atoms:
		indexes = [n for n in mol.graph.neighbors[i]]
		for j in range(0,len(indexes)):
			if not indexes[j] in qm_atoms:
				link_O.append(indexes[j])

'Identify Si atoms attached to link_O in MM region'
link_Si = []
for i in mol.graph.neighbors:
	'loop over all atoms in xyz input file'
	if i in link_O:
		indexes = [n for n in mol.graph.neighbors[i]]
		for j in range(0,len(indexes)):
			if not indexes[j] in qm_atoms and mol.symbols[indexes[j]] == 'Si':
				link_Si.append(indexes[j])

for i in mol.graph.neighbors:
	indexes = [n for n in mol.graph.neighbors[i]]
	if i in qm_atoms:
		type = 0 #no MM information is needed here
	elif i in link_Si:
		type = -24 #different charge is assigned due to H replacement in QM region
		for j in range(0,len(indexes)):
                	if indexes[j] in link_O:
                        	type += 1
	elif mol.symbols[i] == 'O':
		type = -1
	elif mol.symbols[i] == 'Si' or mol.symbols[i] == 'Al':
		type = -2
	elif mol.symbols[i] == 'H':
		type = -3

	'print results'
	if type != 0:
		if i not in link_O:
			x = str(round(mol.coordinates[i][0]/angstrom,3))
			y = str(round(mol.coordinates[i][1]/angstrom,3))
			z = str(round(mol.coordinates[i][2]/angstrom,3))

			print(x+'\t'+y+'\t'+z+'\t'+str(charge_dict[type]))
