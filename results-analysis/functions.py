import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy

'''
Functions complementary to results analysis
'''

def calc_index(index, data_output, exchange, calc_type):
	'''
	finds the entry in data_output and matches it to the item in data_original ['index']
	Input : index in data_original
	exchange: type of exchange functional

	Output: entry in data_output
	'''
	output = 'none' #make default is none, uncless calculation is available
	for item in data_output:
		try:
			if data_output[item]['index'] == index:
				stat1 = calc_type in item
				stat2 = exchange in item
				if stat1 == True:
					if stat2 == True:
						output = item
		except:
			pass

	return output

def Al_Al(atoms):
	'''
	finds distance between 2 Al atoms in qm region
	Input : ase atoms object
	Output: Al-Al distance [unless there is one Al, distance then is 0] and number of Al atoms
	'''

	n_Al = 0 		#number of Al atoms
	Al_index = []	#index of Al atom

	for atom in atoms:
		if atom.symbol == 'Al':
			n_Al += 1
			Al_index.append(atom.index)

	if n_Al == 1:
		distance = 0
	elif n_Al == 2:
		distance = atoms.get_distance(Al_index[0],Al_index[1])
	
	return distance, n_Al

def O_O(atoms, n_Al):
	'''
	finds the distance between two oxygen atoms where H is bonded
	Inputs:
		atoms: ase atoms object
		n_Al : number of Al atoms in qm region
	Output: distance between two oxygen atoms where H is bonded [0 if n_Al = 1]
	'''

	H_index, O_index = [],[]  #indexes of H and O atoms

	if n_Al == 1:
		distance = 0
	elif n_Al == 2:
		for atom in atoms:
			'find H atoms in qm region'
			if atom.symbol == 'H':
				H_index.append(atom.index)

		for H in H_index:
			'find closest O to the H'
			d_min = 100 #start with a large value (dummy value)

			for atom in atoms:
				if atom.symbol == 'O':
					if atoms.get_distance(H, atom.index) < d_min:
						d_min  =  atoms.get_distance(H, atom.index)
						O_atom = atom.index

			O_index.append(O_atom)
	
		distance = atoms.get_distance(O_index[0], O_index[1])
	
	return distance

def Pd_H(atoms, n_Al):
	'''
	finds the distance between Pd and H 
	Inputs:
		atoms: ase atoms object
		n_Al : number of Al atoms in qm region
	Output: distance between H and Pd atom [0 if n_Al = 1]
	'''

	if n_Al == 1:
		distance = 0
	elif n_Al == 2:
		for atom in atoms:
			'find H atoms in qm region'
			if atom.symbol == 'H':
				H_index  = atom.index
			elif atom.symbol == 'Pd':
				Pd_index = atom.index
	
		distance = atoms.get_distance(Pd_index, H_index)
	
	return distance


#def sort(x_pos, E):
'''
	Sort energies of structures from highest to lowest
	Inputs:
		x_pos: label (number) of structures
		E    : energy of the structure
	Outputs:
		new_label: sorted name of each label corresponding to new_E
		new_E    : sorted energies
		x_pts	 : for plotting purposes, from 0 to len(E)
	E_data, new_label, new_E, x_pts = {},[],[],[]
	E_copy = deepcopy(E)

	for index, item in enumerate(E_copy):
		'define a dict with entries being structure name'
		E_data[x_pos[index]] = E_copy[index]
	E_copy.sort() #sort energies from lowest to highest
	for index, E_item in enumerate(E_copy):
		x_pts.append(index)			#x-axis points in the plot
		x = list(E_data.values()).index(E_item) 
		new_label.append(x_pos[x])		#name of the label of each structure
		new_E.append(E_item)			#sorted energies
		
	return new_label, new_E, x_pts
'''
'''
def sort(x_pos, E):
	E_data, new_label, new_E, x_pts = {},[],[],[]
	E_copy = deepcopy(E)

	for index, item in enumerate(E_copy):
		'define a dict with entries being structure name'
		E_data[x_pos[index]] = E_copy[index]
	E_copy.sort() #sort energies from lowest to highest
	for index, E_item in enumerate(E_copy):
		x_pts.append(index)			#x-axis points in the plot
		x = list(E_data.values()).index(E_item) 
		new_label.append(x_pos[x])		#name of the label of each structure
		new_E.append(E_item)			#sorted energies
		
	return new_label, new_E, x_pts
'''


def sort(x_pos, E):
	'''
	Sort energies of structures from highest to lowest
	Inputs:
		x_pos: label (number) of structures
		E    : energy of the structure
	Outputs:
		new_label: sorted name of each label corresponding to new_E
		new_E    : sorted energies
		x_pts	 : for plotting purposes, from 0 to len(E)
	'''
	E_data, new_label, new_E, x_pts = {},[],[],[]
	E_copy = deepcopy(E)
	'''
	for index, item in enumerate(E_copy):
		'define a dict with entries being structure name'
		E_data[x_pos[index]] = E_copy[index]
	E_copy.sort() #sort energies from lowest to highest
	for index, E_item in enumerate(E_copy):
		x_pts.append(index)			#x-axis points in the plot
		x = list(E_data.values()).index(E_item)
		new_label.append(x_pos[x])		#name of the label of each structure
		new_E.append(E_item)			#sorted energies
	'''

	for index, item in enumerate(E_copy):
		'define a dict with entries being structure name'
		E_data[x_pos[index]] = E_copy[index]
	E_copy.sort() #sort energies from lowest to highest
	for index, E_item in enumerate(E_copy):
		x_pts.append(index)			#x-axis points in the plot
		x = list(E_data.values()).index(E_item)
		new_label.append(x_pos[x])		#name of the label of each structure
		new_E.append(E_item - min(E_copy))			#sorted energies

	return new_label, new_E, x_pts


def min_H(ref, H_data,calc_type):
	'''
	out of the 16 possible configurations of H sites, identify the lowest energy
	Inputs: ref - name of the original zeolite from which the Pd2+ was created (as well as H2+ by definition)
	H_data: dir where compensating protons are saved
	Output: minimum energy of the zeolite structure with protons (out of the 16 possibilities)
	'''

	'load data of H adsorbed on zeolite'
	with open(H_data+'/data.json','r') as read_file:
		data_H = json.load(read_file)

	with open(H_data+'/data_output.json','r') as read_file:
		data_H_output = json.load(read_file)

	list_ref = []	#save items that share the same reference

	for item in data_H:
		'generate a list of items sharing the same zeolite reference (from which H2+ is created)'
		if data_H[item]['reference'] == ref:
			if item != ref:
				list_ref.append(item[0:-5])

	min_energy = 0	#define a high energy as starting point
	for item in data_H_output:
		'find energy of each item in list of references'
		if 'traj' in item:
			'avoids wrong entries such as incomplete or tmp'
			if data_H_output[item]['index'] in list_ref:
				if calc_type in item:
					'only extract calc from sp'
					if data_H_output[item]['energy'] < min_energy:
						a = item
						min_energy =  data_H_output[item]['energy']

	return min_energy

def rxn_energy(E, zeolite_H, oxidation):
	'''
	calculates rxn energy based on the following rxn:
	Pd + [2H(+) z(2-)] --> Pd(2+) Z(2-) + H2O - 1/2 O2
	oxidation: oxidation state of Pd [1 or 2]
	'''
	Pd  = -127.913481461
	H2O = -76.439413334
	O2  = -150.2765115625

	if oxidation == 2:
		energy = E + H2O - 0.5*O2 - Pd - zeolite_H
	elif oxidation == 1:
		energy = E + 0.5*H2O - 0.25*O2 - Pd - zeolite_H

	return energy*27.2114
