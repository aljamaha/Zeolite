import json, os
import matplotlib.pyplot as plt
from ase import atoms, io
from copy import deepcopy

'''
Results Analysis
'''

'Inputs'
sorted_plot = True		#if True, bar plots of energies is sorted from lowest to highest
#ring_color  = False		#if True, bar plots with ring type
data_dir    = '/home/aljama/CHA-Pd/data/'			#dir where json data are saved
calc_dir    = '/home/aljama/CHA-Pd/calculations/'		#dir where calculations are done
results_dir = '/home/aljama/CHA-Pd/results-analysis/' 	#dir where results are to be saved
H_data      = '/home/aljama/CHA/data/'	#dir where data for H adsorptions sites are saved
ring_data   = '/home/aljama/CHA/data/' 	#dir with information on ring types, etc.

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

def sort(x_pos, E, colour):
	'''
	Sort energies of structures from highest to lowest
	Inputs:
		x_pos : label (number) of structures

		E     : energy of the structure
		colour: list of colors for each data based on ring type
	Outputs:
		new_label: sorted name of each label corresponding to new_E
		new_E    : sorted energies
		x_pts	 : for plotting purposes, from 0 to len(E)
		new_color: sorted list of colors
	'''
	E_data, new_label, new_E, x_pts,new_color = {},[],[],[],[]
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
		new_color.append(colour[x])
		
	return new_label, new_E, x_pts, new_color

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

def find_MR(ref, ring_data, color_dict, original):
	'''
	finds the type of MR based on the reference
	Inputs:
		ref      : the original reference for the calculation
		ring_data: dict of each type of MR and the original zeolites it represents
		color_dict: dict of each type of MR and the color associated with eac
		original: index name of the calculations conducted (under calculations folder)
		
	Output: color of the index calculation based on the ref
	'''
	match_index = ''
	index = ref[0:-5] #removes traj part
	for item in ring_data:
		for i in ring_data[item]:
			if str(i) == str(index):
				match_index = item
				break
	if match_index != '':
		try:
			output = color_dict[match_index]
		except:
			output = 'b'
	else:
		output = 'y'

	return output

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
label, E, Al_d, x_pos, colour, E_min, x_min, label_min = {}, [], [],[],[],[],[],[]	 #saving values for all entries
Al_d, E_d, c_dict = {},{},{}
first_item = True

'add color input based on the type of MR'
try:
		with open(ring_data+'/CHA_ring_type.json', 'r') as read:
			ring_data = json.load(read)		#load data
		color = ['b','r','g','c','k']		#assign color names
		#color = ['b','b','b','b','b','b']		#assign color names

		color_dict = {}				
		for index, item in enumerate(ring_data):
			color_dict[item] = color[index]		#assign color to each MR
except:	
		color_dict = {}
		print('Missing json data file for ring data')

'manual deletions of repeated refernces'
print('Warning! Manual deletion of some references!')
del references['6.traj'] 	#repeated 4 MR
del references['5.traj']	#repeated as 3
del references['11.traj']	#repeated as 17 (stackd 6 MR)
del references['8.traj']	#repeated as 7
del references['9.traj']	#repeated as 7
del references['10.traj']	#repeated as 7
del references['23.traj']	#repeated as 25
del references['26.traj']	#repeated as 25

for ref in references:
	'loop over each reference'

	if os.path.isdir(results_dir+'/traj') == False:
		'create folder if it does not exist'
		os.system('mkdir '+results_dir+'/traj')

	if os.path.isdir(results_dir+'/traj-surroundings') == False:
		'create folder if it does not exist'
		os.system('mkdir '+results_dir+'/traj-surroundings')

	E_tmp, x_tmp, label_tmp = [], [], []
	
	for item in references[ref]:
		'each item under reference'

		index = item[0:-5] 			#remves .traj from the name
		data_output_entry = calc_index(index) 	#check corresponding name in data_output
		if data_output_entry != 'none':		#check calcuation dir is available
			if data_output[data_output_entry]['status'] == 'complete':
				'check calc is completed, then copy traj files to new folder'
				os.system('cp '+calc_dir+'/'+data_output_entry+'/qm-final.xyz '+results_dir+'/traj/'+item[0:-5]+'.traj')
				os.system('cp '+calc_dir+'/'+data_output_entry+'/surroundings-final.traj '+results_dir+'/traj-surroundings/'+item[0:-5]+'.traj')

				'calculate rxn energy'
				ref_H = data_original[item]['reference'] #reference for 16 H calculations
				zeolite_H = min_H(ref_H)
				E_qmmm = data_output[data_output_entry]['energy']
				E_rxn = rxn_energy(E_qmmm, zeolite_H)
				if first_item == True:
					E_ref = E_rxn
				E.append(E_rxn - E_ref)
				first_item = False

				'Al-Al distance'
				atoms = io.read(calc_dir+data_output_entry+'/qm-initial.traj')
				Al_distance, n_Al = Al_Al(atoms)
				data_output[data_output_entry]['Al-Al distance'] = round(Al_distance,3)
				Al_d[index], E_d[index], label[index] =  Al_distance, E[-1], index
				c_dict[index] = find_MR(ref, ring_data, color_dict,index)
				
				'Accumulate data'
				x_pos.append(int(index))	#x-asis position
				E_tmp.append(E[-1])
				x_tmp.append(int(index))
				label_tmp.append(index)
				c = index			#to avoid confusion since index is called again later

	try:
		'finding only the minimum out of each reference'
		minimum_Energy = min(E_tmp)

		for index, item in enumerate(E_tmp):
			if item == minimum_Energy:
				E_min.append(item)
				x_min.append(x_tmp[index])
				label_min.append(label_tmp[index])
				break	

		'color code for bar plots'
		try:
			colour.append(find_MR(ref, ring_data, color_dict, c))
		except:
			colour.append('b')	

	except:
		print('reference '+ref+' does not contain values')

print('Ring data :', ring_data)
print('Color dict:', color_dict)

'''Plots'''
if sorted_plot == True:
	'bar plot (sorted)'
	print('WARNING: manual entry of NNN/NNNN')
	new_x, new_E, x_pts, new_C = sort(x_min, E_min, colour)
	for index, item in enumerate(new_E):
		#print(x_min[index], E_min[index])
		plt.bar(x_pts[index], new_E[index],color=new_C[index], align='center', alpha=1)
		if new_x[index]>45:
			plt.text(x_pts[index]-0.25,0.5,'NNNN', rotation = 90)
		else:
			plt.text(x_pts[index]-0.25,0.5,'NNN', rotation = 90)
	plt.xticks(x_pts, new_x, rotation = 90)
	plt.ylabel('Energy (eV)')
	plt.show()
else:
	'bar plot (not sorted)'
	for index, item in enumerate(E_min):
		plt.bar(x_min[index], E_min[index], color=colour[index], align='center', alpha=1)
	plt.xticks(x_min, x_min, rotation=90)
	plt.ylabel('Energy (eV)')
	plt.show()

'Al-Distance plot'
for item in new_x:
	try:
		plt.plot(Al_d[str(item)], E_d[str(item)], c_dict[str(item)]+'o', markersize=6)
	except:
		plt.plot(Al_d[str(item)], E_d[str(item)], 'bo', markersize=6)
	plt.text(Al_d[str(item)], E_d[str(item)], label[str(item)])
plt.xlabel('Al-Al Distance (A)', fontsize = 10)
plt.ylabel('Energy (eV)', fontsize = 10)
plt.show()

'save data with Al-Al distance'
with open(data_dir+"data_output.json", "w") as write_file:
    json.dump(data_output, write_file, indent=4)
