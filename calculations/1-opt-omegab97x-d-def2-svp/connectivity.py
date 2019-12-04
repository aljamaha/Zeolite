#import packages
import sys
from molmod import *	#further info are here: http://molmod.github.io/molmod/tutorial/install.html

'''
Objective: convert an xyz coordinates to qmmm molecule format
Input format: python qm-mm-connectivity <xyz coordinates> <List of qm atoms>
Output: $molecule section in qm-mm
'''

'Modified version from Jeroen Van Der Mynsbrugge original script'
'''
if len(sys.argv) !=3:
	'raise an error when less than two arguments are provided'
	sys.stderr.write("Wrong input format\nInput format: python qm-mm-connectivity.py <xyz file> <#of qm atoms>\n")
	quit()
'''
'qm_atoms input is a list. Convert it into a numeric python list'
qm_atoms = list(sys.argv[2:])	#List number of qm atoms
empty    = 0			#stores the number of empty entries into the list

for index, item in enumerate(qm_atoms):
	if '[' in item:	
		item = item.replace("[",'')
	if ']' in item:
		item = item.replace("]",'')
	if ',' in item:
		item = item.replace(",",'') 
	if item == "":
		empty +=1 
		qm_atoms[index] = item
	else:
		qm_atoms[index] = int(item)

for i in range(0,empty):
	qm_atoms.remove("")

mol = Molecule.from_file(sys.argv[1])	#reads the input <xyz coordinates>
mol.set_default_masses()		#reads the mass of each atom in xyz input
assert(mol.graph is None)
mol.set_default_graph()			#derive a molecular graph based on geometry

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
	print("%s\t % f\t % f\t % f\t % i\t %i\t %i\t %i\t %i\t" % (str(mol.symbols[i]), mol.coordinates[i][0]/angstrom, mol.coordinates[i][1]/angstrom, mol.coordinates[i][2]/angstrom, type, indexes[0]+1 if len(indexes) > 0 else 0, indexes[1]+1 if len(indexes) > 1 else 0, indexes[2]+1 if len(indexes) > 2 else 0, indexes[3]+1 if len(indexes) > 3 else 0))


