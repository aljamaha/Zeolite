#import packages
import sys, json
from molmod import *	#further info are here: http://molmod.github.io/molmod/tutorial/install.html
from ase import io

'''
Objective: convert an xyz coordinates to QM/MM molecule format
Input format: python qm-mm-connectivity <xyz coordinates> <List of qm atoms> <Zeolite name (OPTIONAL)>
Output: $molecule section in QM/MM calculations
'''

'Modified version from Jeroen Van Der Mynsbrugge original script'

try:
	'Input arguments. qm_atom is a list converted into a numeric python list'
	mol = Molecule.from_file(sys.argv[1])	#reads the input <xyz coordinates>
	qm_atoms  = list(sys.argv[2:-1])		#List of qm atoms
except:
	print('Wrong input format!\npython connectivity.py <xyz-file> <list of qm atoms> <Zeolite name (OPTIONAL)')
	exit()

try:
	atoms = io.read(sys.argv[1][0:-4]+'.traj')
except:
	print('input traj file not available')
	exit()

try:
	zeolite  = sys.argv[-1]	
	with open(zeolite+"-NL.json", "r") as read_file:
		data = json.load(read_file)
	n_atoms = len(data)
	input_Zeolite = True
	#qm_atoms = qm[0:-1]
except:
	input_Zeolite = False
	#qm_atoms = qm

def clean_qm_atoms(qm_atoms):
	'converts list of qm_atoms to a str easily used later on'
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

	return qm_atoms

qm_atoms = clean_qm_atoms(qm_atoms)
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
	if input_Zeolite == False:
		print("%s\t % f\t % f\t % f\t % i\t %i\t %i\t %i\t %i\t" % (str(mol.symbols[i]), mol.coordinates[i][0]/angstrom, mol.coordinates[i][1]/angstrom, mol.coordinates[i][2]/angstrom, type, indexes[0]+1 if len(indexes) > 0 else 0, indexes[1]+1 if len(indexes) > 1 else 0, indexes[2]+1 if len(indexes) > 2 else 0, indexes[3]+1 if len(indexes) > 3 else 0))
	elif input_Zeolite == True:
		if i < n_atoms:
			print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( str(mol.symbols[i]), round(mol.coordinates[i][0]/angstrom,3), round(mol.coordinates[i][1]/angstrom,3), round(mol.coordinates[i][2]/angstrom,3), type, data[i][1], data[i][2], data[i][3], data[i][4])  )
		else:
			print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( str(mol.symbols[i]), round(mol.coordinates[i][0]/angstrom,3), round(mol.coordinates[i][1]/angstrom,3), round(mol.coordinates[i][2]/angstrom,3), type, 0, 0, 0, 0)  )




