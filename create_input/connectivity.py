#import packages
import sys, json
from molmod import *	#further info are here: http://molmod.github.io/molmod/tutorial/install.html
from ase import io

'Modified version from Jeroen Van Der Mynsbrugge original script'

def connect(xyz_file, qm_atoms, output_name = 'tmp', zeolite=None):
	'''
	Objective: convert an xyz coordinates to QM/MM molecule format
	Inputs:
		xyz_file    : name of the xyz file
		qm_atoms    : list of atoms in the qm region
		zeolite	    : name of the zeolite atoms
		output_name : name of the output file (default is tmp)
	Output: $molecule section in QM/MM calculations
	'''


	f = open(output_name, 'w')
	mol       = Molecule.from_file(xyz_file)	#reads the input <xyz coordinates>
	atoms     = io.read(xyz_file)
	
	if zeolite != None:
		with open(zeolite+"-NL.json", "r") as read_file:
			data = json.load(read_file)
		n_atoms = len(data)

	#qm_atoms = clean_qm_atoms(qm_atoms)
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
		if zeolite == None:
			f.write("%s\t % f\t % f\t % f\t % i\t %i\t %i\t %i\t %i\t" % (str(mol.symbols[i]), mol.coordinates[i][0]/angstrom, mol.coordinates[i][1]/angstrom, mol.coordinates[i][2]/angstrom, type, indexes[0]+1 if len(indexes) > 0 else 0, indexes[1]+1 if len(indexes) > 1 else 0, indexes[2]+1 if len(indexes) > 2 else 0, indexes[3]+1 if len(indexes) > 3 else 0)+'\n')
		else: 
			if i < n_atoms:
				f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( str(mol.symbols[i]), round(mol.coordinates[i][0]/angstrom,3), round(mol.coordinates[i][1]/angstrom,3), round(mol.coordinates[i][2]/angstrom,3), type, data[i][1], data[i][2], data[i][3], data[i][4])+'\n')
			else:
				f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format( str(mol.symbols[i]), round(mol.coordinates[i][0]/angstrom,3), round(mol.coordinates[i][1]/angstrom,3), round(mol.coordinates[i][2]/angstrom,3), type, 0, 0, 0, 0)+'\n'  )
	f.close()

