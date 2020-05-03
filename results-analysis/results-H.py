import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy
import numpy as np

'''
Results Analysis (energy and O-O distances)
add a description here
'''

'Inputs'
plotting    	    = True		#if True, plot results for each reference structure
plot_individual	    = True		#show individual plots for each reference (as a function of O-O distance)
sorted_plot	    = False		#if True, bar plots of energies is sorted from lowest to highest
not_sorted	    = False		#make a plot of all data aggregated (O-distance vs. Energy)
plotting_overall    = True		#if True, make an overall plot of all results
label_pts	    = False		#add labels to pts in scatter plt
Al_shading	    = True		#add vertical line referring to Al-Al distance
data_dir    = '/home/aljama/CHA-full-MR/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/CHA-full-MR/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/CHA-full-MR/results-analysis/' 	#dire where results are to be saved
#colors      = ['g','b','r','k','y','c','m']

'Load data from json files'
with open(data_dir+"data_output.json", "r") as read_file:
    data_output = json.load(read_file)

with open(data_dir+"data.json", "r") as read_file:
    data_original = json.load(read_file)

def calc_index(index):
	'''
	finds the entry in data_output and matches it to the item in data_original ['index']
	Input : index in data_original
	Output: entry in data_output
	'''
	output = 'none' #make default is none, uncless calculation is available

	for item in data_output:
		if data_output[item]['index'] == index:
			stat = '-sp-' in item
			if stat == True:
				output = item
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


def min_H(ref):
	'''
	out of the 16 possible configurations of H sites, identify the lowest energy
	Inputs: ref - name of the original zeolite from which the Pd2+ was created (as well as H2+ by definition)
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
		if data_H_output[item]['index'] in list_ref:
			if '-sp-' in item:
				'only extract calc from sp'
				if data_H_output[item]['energy'] < min_energy:
					a = item
					min_energy =  data_H_output[item]['energy']

	return min_energy

def rxn_energy(E, zeolite_H):
	'''
	calculates rxn energy based on the following rxn:
	Pd + [2H(+) z(2-)] --> Pd(2+) Z(2-) + H2O - 1/2 O2
	'''
	Pd  = -1496.0850202371
	H2O = -76.439413334
	O2  = -150.2765115625

	energy = E + H2O - 0.5*O2 - Pd - zeolite_H

	return energy*27.2114

'accumulate reference entries (templates from which calculations were created and run)'
references = {} #references for data

for item in data_original:
	'entries in references dictionary'
	if data_original[item]['reference'] not in references:
		references[data_original[item]['reference']] = []

for ref in references:
	'name of folders that share same reference'
	for item in data_original:
		if data_original[item]['reference'] == ref:
			if ref != item:
				references[ref].append(item)

'accumulate traj files'
O_d_all, E_all = [],[] #those are for all calculations to get a comprehensive picture

for ref in references:
	'loop over each reference'
	x_pos, E, first_item, O_d, label = [],[], True, [],[]	#For plotting purposes

	if os.path.isdir(results_dir+ref) == False:
		'create folder if it does not exist'
		os.system('mkdir '+results_dir+ref)

	os.chdir(results_dir+ref)

	for item in references[ref]:
		'each item under reference'
		index = item[0:-5] #remves .traj from the name
		data_output_entry = calc_index(index) #check corresponding name in data_output
		if data_output_entry != 'none':
			'check calcuation dir is available'
			if data_output[data_output_entry]['status'] == 'complete':
				'check calc is completed, then copy traj files to new folder'
				os.system('cp '+calc_dir+'/'+data_output_entry+'/qm-initial.traj '+item[0:-5]+'.traj')
				x_pos.append(int(index))	#x-asis position
				if first_item == True:
					E_ref = data_output[data_output_entry]['energy']
				E.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) #convert from Hartree to e. values of energies for y-axis
				E_all.append( (data_output[data_output_entry]['energy']- E_ref)*27.2114 ) 
				first_item = False

				'Al-Al distance'
				atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
				Al_distance, n_Al = Al_Al(atoms)
				data_output[data_output_entry]['Al-Al distance'] = round(Al_distance,3)

				'O-O distance'
				O_O_distance = O_O(atoms, n_Al)
				data_output[data_output_entry]['O-O distance'] = round(O_O_distance,3)
				O_d.append( round(O_O_distance,3) )
				O_d_all.append( round(O_O_distance,3) )
				
				label.append(index)

	if plotting == True:	
		if sorted_plot == True:
			'bar plot (sorted)'
			new_x, new_E, x_pts = sort(x_pos, E)
			plt.bar(x_pts, new_E, align='center', alpha=1)
			plt.xticks(x_pts, new_x)
			plt.ylabel('Energy (eV)')
			plt.show()

		elif not_sorted == True:
			'bar plot (not sorted)'
			plt.bar(x_pos, E, align='center', alpha=1)
			plt.xticks(x_pos, x_pos)
			plt.ylabel('Energy (eV)')
			plt.show()
		if plotting_overall == True:
			try:
				tmp, tmp_ind = np.min(E), E.index(min(E))
				plt.plot(O_d[tmp_ind], tmp - np.min(E) , 'sk', markersize=10)
				plt.plot(O_d, E - np.min(E), 'o', markersize=6,label=ref[0:-5])
				plt.xlabel('O-O Distance (A)', fontsize = 14)
				plt.ylabel('Energy (eV)', fontsize = 14)
				plt.tick_params(labelsize=14)
				if Al_shading == True:
					plt.plot([Al_distance,Al_distance],[-2,2],'r--')
					plt.plot([Al_distance-0.8,Al_distance-0.8],[-2,2],'k--')
					plt.plot([Al_distance+0.8,Al_distance+0.8],[-2,2],'k--')
				if label_pts == True:
					for i, lab in enumerate(label):
						plt.text(O_d[i], E[i], lab)
				if plot_individual == True:
					plt.show()
			except:
				pass

if plotting_overall == True:
	plt.legend()
	plt.xlim([1,10])
	plt.show()
